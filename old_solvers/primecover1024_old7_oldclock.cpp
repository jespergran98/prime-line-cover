/*
 * primecover1024.cpp
 *
 * Exact minimum line cover solver for prime points (i, p_i).
 *
 * Core model:
 *   - Heavy sets are all collinear lines containing at least 3 points.
 *   - The DFS chooses a minimum-cardinality subset of heavy lines.
 *   - Any residual uncovered points are closed exactly with 2-point / 1-point
 *     lines, costing ceil(residual / 2).
 *
 * Search strategy:
 *   - 1024-ready bit-parallel masks for all coverage operations.
 *   - Binary include/exclude branching on one gain-positive heavy line at a
 *     time, ordered by current gain and local overlap pressure.
 *   - Incremental state tracking: productive_count / size_freq updated via
 *     delta application, with per-depth state-frame reuse to avoid repeated
 *     heap allocation in the monolithic DFS.
 *   - Parallel exact frontier splitting for monolithic instances: worker
 *     threads explore disjoint DFS subtrees while sharing one atomic
 *     incumbent so any improvement prunes every other thread immediately.
 *   - Optimistic lower bound + branch-and-bound pruning.
 *   - Memoization on (active mask, dangling parity) in the monolithic solver.
 *   - Exact component fallback when the active heavy-line graph truly
 *     decomposes; dense rigid instances stay on the faster monolithic path.
 *
 * Experiment log:
 *   1. Replaced the active ExactGainSolver with the older direct-cost
 *      monolithic branch-and-bound family, adding exact child pre-pruning and
 *      duplicate-child elimination to make that swap as competitive as
 *      possible. The result was still dramatically worse on this machine: the
 *      required 1..270 sweep stopped at N=181 after N=181 took 50.931s, with a
 *      1..181 total of 295.248714s versus the gain solver's full 1..270 totals
 *      of 59.250765s and 60.120906s. The simplest explanation is that the
 *      direct cost search still explores far more states before it can certify
 *      infeasibility on these dense monolithic instances, so future work
 *      should not replace the gain-based engine with this solver family
 *      without a materially stronger root-to-leaf proof bound.
 *   2. Replaced unconstrained exact-gain maximization with an incumbent-beating
 *      decision search that carried an explicit heavy-line budget and used the
 *      budgeted LP dual for additional gain as its pruning bound. It stayed
 *      correct but regressed sharply once the hard tail began: compared with
 *      the baseline, N=246 went from 0.506134s to 7.375714s, N=247 from
 *      0.574438s to 8.416809s, and N=248 from 0.382636s to well over 8s before
 *      the run was aborted. The simplest explanation is that proving
 *      infeasibility separately for successive target costs re-explores too
 *      much of the same availability-state search tree, so future work should
 *      not switch to target-by-target decision passes unless it also brings a
 *      much stronger cross-target memoisation scheme or a materially smaller
 *      decision tree.
 *   3. Kept the original one-pass exact-gain DFS but added incumbent-aware
 *      pruning from the same budgeted LP dual, so every node was bounded by
 *      the number of heavy lines still available to beat the current best
 *      cost. That avoided the repeated target passes from experiment 2, but it
 *      still regressed badly in the hard tail: N=244 rose from 0.122117s to
 *      1.736302s, N=245 from 0.312369s to 1.261494s, and N=246 from 0.506134s
 *      to 7.195130s before the run was aborted. The simplest explanation is
 *      that the current availability-state partition tree does not benefit
 *      enough from a line-budget dual to repay the extra bound cost, so future
 *      work should not graft budget-aware LP pruning onto this DFS without a
 *      cheaper relaxation or a different state decomposition.
 *   4. Added a direct remaining-cost lower bound to the live ExactGainSolver
 *      by adapting the legacy size-histogram and overlap-dual cost bounds and
 *      pruning any branch whose chosen lines plus that bound could not beat
 *      the incumbent cost. The root measurements showed the new bound was
 *      actually weaker than the existing gain-side Lagrangian bound on the
 *      hard tail (for example N=274 stayed at cost_lb=40 from the live bound
 *      but only cost_lb_direct=36 from the new one), and the full run regressed
 *      enough to stop at N=277 after 15.467s. The simplest explanation is that
 *      this cost relaxation throws away too much of the gain model's shared
 *      coverage structure while still charging extra per-node work, so future
 *      attempts should not front-load the gain solver with the legacy cost
 *      bounds unless they first prove dominance over the existing root bound.
 *   5. Retried the same direct-cost bound integration after fixing the
 *      accidental pair-closure prune from the first attempt, so the added
 *      histogram/overlap-dual bound was now sound and the solver again matched
 *      the oracle through N=290. The hard-tail root bound was still weaker
 *      than the live gain proof at every tested root (N=274..290 stayed at
 *      direct_cost_lb=36..38 versus gain_cost_lb=40..44), and the extra bound
 *      work still regressed the full tail badly: in the diagnostic run N=284
 *      rose from 3.387437s to 9.975613s and N=290 from 5.866359s to 16.849s.
 *      The simplest explanation is that this relaxed cost model loses exactly
 *      the overlap structure the gain bound already captures, so future work
 *      should not spend more node time on the legacy cost dual unless a new
 *      relaxation first proves tighter at the root on the hard tail.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <queue>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int    kBitCapacity          = 1024;
constexpr int    kBitWords             = kBitCapacity / 64;
constexpr int    kStartN               = 650;   // N to start from in the timed sweep; set to 0 to run from 1
constexpr int    kExecutionLimit       = 720;   // hard cap and default run-to N
// Maximum wall-clock seconds allowed per N. Set to 0.0 to disable (no limit). 60=1 min
constexpr double kPerNTimeLimitSeconds = 15;

// ---------------------------------------------------------------------------
// BitMask1024
// ---------------------------------------------------------------------------
struct BitMask1024 {
    std::array<std::uint64_t, kBitWords> words{};

    void set(int bit) noexcept   { words[bit >> 6] |=  (1ULL << (bit & 63)); }
    void reset(int bit) noexcept { words[bit >> 6] &= ~(1ULL << (bit & 63)); }
    bool test(int bit) const noexcept {
        return ((words[bit >> 6] >> (bit & 63)) & 1ULL) != 0;
    }

    int count(int used_words) const noexcept {
        int total = 0;
        for (int i = 0; i < used_words; ++i) total += std::popcount(words[i]);
        return total;
    }

    bool intersects(const BitMask1024& other, int used_words) const noexcept {
        for (int i = 0; i < used_words; ++i)
            if ((words[i] & other.words[i]) != 0) return true;
        return false;
    }

    int intersection_count(const BitMask1024& other, int used_words) const noexcept {
        int total = 0;
        for (int i = 0; i < used_words; ++i)
            total += std::popcount(words[i] & other.words[i]);
        return total;
    }

    BitMask1024 and_not(const BitMask1024& other) const noexcept {
        BitMask1024 out;
        for (int i = 0; i < kBitWords; ++i) out.words[i] = words[i] & ~other.words[i];
        return out;
    }

    static BitMask1024 prefix(int n) noexcept {
        BitMask1024 out;
        const int full_words = n / 64;
        for (int i = 0; i < full_words; ++i)
            out.words[i] = std::numeric_limits<std::uint64_t>::max();
        if (const int tail = n % 64; tail != 0 && full_words < kBitWords)
            out.words[full_words] = (1ULL << tail) - 1ULL;
        return out;
    }

    template <class Fn>
    void for_each_set_bit(int used_words, Fn&& fn) const {
        for (int wi = 0; wi < used_words; ++wi) {
            std::uint64_t word = words[wi];
            while (word != 0) {
                fn((wi << 6) + std::countr_zero(word));
                word &= word - 1;
            }
        }
    }

    friend BitMask1024 operator&(const BitMask1024& a, const BitMask1024& b) noexcept {
        BitMask1024 out;
        for (int i = 0; i < kBitWords; ++i) out.words[i] = a.words[i] & b.words[i];
        return out;
    }
};

// ---------------------------------------------------------------------------
// Keys and hashes
// ---------------------------------------------------------------------------
struct LineKey {
    std::int64_t a{}, b{}, c{};
    bool operator==(const LineKey&) const = default;
};

struct LineKeyHash {
    static std::uint64_t mix(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    std::size_t operator()(const LineKey& k) const noexcept {
        std::uint64_t h = mix(static_cast<std::uint64_t>(k.a));
        h ^= mix(static_cast<std::uint64_t>(k.b) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        h ^= mix(static_cast<std::uint64_t>(k.c) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        return static_cast<std::size_t>(h);
    }
};

// Used by ComponentProfileSolver (memoises without dangling parity).
struct StateKey {
    std::array<std::uint64_t, kBitWords> words{};
    bool operator==(const StateKey&) const = default;
};

struct StateKeyHash {
    static inline int used_words = kBitWords;

    static std::uint64_t mix(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    std::size_t operator()(const StateKey& k) const noexcept {
        std::uint64_t h = mix(0x123456789abcdef0ULL);
        for (int i = 0; i < used_words; ++i)
            h ^= mix(k.words[static_cast<std::size_t>(i)] +
                     0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        return static_cast<std::size_t>(h);
    }
};

// Used by MonolithicExactSolver (memoises with dangling parity).
struct StateKeyWithParity {
    std::array<std::uint64_t, kBitWords> words{};
    bool dangling = false;
    bool operator==(const StateKeyWithParity&) const = default;
};

struct StateKeyWithParityHash {
    static inline int used_words = kBitWords;

    static std::uint64_t mix(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
    std::size_t operator()(const StateKeyWithParity& k) const noexcept {
        std::uint64_t h = mix(static_cast<std::uint64_t>(k.dangling));
        for (int i = 0; i < used_words; ++i)
            h ^= mix(k.words[static_cast<std::size_t>(i)] +
                     0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
        return static_cast<std::size_t>(h);
    }
};

// ---------------------------------------------------------------------------
// Core data structures
// ---------------------------------------------------------------------------
struct HeavyLine {
    BitMask1024      mask;
    std::vector<int> points;
    int              activate_at = 0;
};

// Incremental state threaded through the monolithic DFS.
// All fields are updated via apply_delta, and per-depth buffers are reused
// so sibling branches avoid repeated heap allocation.
struct IncrementalState {
    std::vector<std::uint16_t> line_cover;       // current active intersection size per line
    std::vector<std::uint16_t> productive_count; // # productive lines covering each point
    std::vector<int>           size_freq;        // size_freq[k] = # productive lines with cover == k
};

struct ActiveComponent {
    std::vector<int> line_ids;
    BitMask1024      point_mask;
    int              point_count = 0;
};

using Summary2 = std::array<int, 2>;   // {even-parity cost, odd-parity cost}

struct GreedyCoverSolution {
    int              cost = 0;
    std::vector<int> line_ids;
};

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------
[[nodiscard]] int words_for_n(int n) noexcept { return (n + 63) / 64; }
[[nodiscard]] int ceil_half(int v)   noexcept { return (v + 1) / 2; }

[[nodiscard]] std::string format_seconds(double s) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(s < 10.0 ? 6 : 3) << s << 's';
    return out.str();
}

[[nodiscard]] int capped_add(int lhs, int rhs, int cap) noexcept {
    if (lhs > cap || rhs > cap || lhs > cap - rhs) return cap + 1;
    return lhs + rhs;
}

// ---------------------------------------------------------------------------
// Prime generation (uses prime number theorem to avoid undershooting the sieve)
// ---------------------------------------------------------------------------
[[nodiscard]] std::vector<int> generate_primes(int count) {
    if (count <= 0) return {};
    const double n = static_cast<double>(count);
    const int limit = (count < 10)
        ? 64
        : static_cast<int>(n * (std::log(n) + std::log(std::log(n)) + 3.0)) + 256;
    std::vector<bool> is_p(static_cast<std::size_t>(limit + 1), true);
    is_p[0] = is_p[1] = false;
    for (int p = 2; p * p <= limit; ++p)
        if (is_p[static_cast<std::size_t>(p)])
            for (int q = p * p; q <= limit; q += p)
                is_p[static_cast<std::size_t>(q)] = false;
    std::vector<int> primes;
    primes.reserve(static_cast<std::size_t>(count));
    for (int v = 2; v <= limit && static_cast<int>(primes.size()) < count; ++v)
        if (is_p[static_cast<std::size_t>(v)]) primes.push_back(v);
    return primes;
}

[[nodiscard]] std::int64_t abs64(std::int64_t x) noexcept { return x < 0 ? -x : x; }

[[nodiscard]] LineKey canonical_line(int x1, int y1, int x2, int y2) {
    std::int64_t a = static_cast<std::int64_t>(y1) - static_cast<std::int64_t>(y2);
    std::int64_t b = static_cast<std::int64_t>(x2) - static_cast<std::int64_t>(x1);
    std::int64_t c = -(a * static_cast<std::int64_t>(x1) + b * static_cast<std::int64_t>(y1));
    const std::int64_t g = std::gcd(std::gcd(abs64(a), abs64(b)), abs64(c));
    if (g) { a /= g; b /= g; c /= g; }
    if (a < 0 || (a == 0 && b < 0) || (a == 0 && b == 0 && c < 0)) {
        a = -a; b = -b; c = -c;
    }
    return {a, b, c};
}

[[nodiscard]] std::vector<HeavyLine> enumerate_heavy_lines(const std::vector<int>& primes) {
    const int n          = static_cast<int>(primes.size());
    const int used_words = words_for_n(n);

    std::unordered_map<LineKey, BitMask1024, LineKeyHash> line_masks;
    line_masks.reserve(static_cast<std::size_t>(n * n));
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            BitMask1024& m = line_masks[canonical_line(i, primes[i], j, primes[j])];
            m.set(i); m.set(j);
        }

    std::vector<HeavyLine> lines;
    lines.reserve(line_masks.size());
    for (auto& [key, mask] : line_masks) {
        if (mask.count(used_words) < 3) continue;
        HeavyLine hl;
        hl.mask = mask;
        mask.for_each_set_bit(used_words, [&](int p) { hl.points.push_back(p); });
        hl.activate_at = hl.points[2] + 1;
        lines.push_back(std::move(hl));
    }
    // Deterministic ordering: first point asc, then size desc, then lexicographic.
    std::sort(lines.begin(), lines.end(), [](const HeavyLine& a, const HeavyLine& b) {
        if (a.points.front() != b.points.front()) return a.points.front() < b.points.front();
        if (a.points.size()  != b.points.size())  return a.points.size()  > b.points.size();
        return a.points < b.points;
    });
    return lines;
}

[[nodiscard]] GreedyCoverSolution greedy_upper_bound_solution(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
{
    const int words = words_for_n(n);
    BitMask1024 active = BitMask1024::prefix(n);
    GreedyCoverSolution result;
    int count = n;
    while (true) {
        int best = -1, bc = 2;
        for (int id : active_ids) {
            const int c = lines[id].mask.intersection_count(active, words);
            if (c > bc || (c == bc && best != -1 &&
                           lines[id].points.size() > lines[best].points.size()))
                { bc = c; best = id; }
        }
        if (best < 0) break;
        active = active.and_not(lines[best].mask);
        count -= bc;
        ++result.cost;
        result.line_ids.push_back(best);
    }
    result.cost += ceil_half(count);
    return result;
}

struct GainSolveResult {
    int              best_gain = 0;
    std::vector<int> chosen_line_ids;
};

[[nodiscard]] int cost_from_gain(int n, int gain) noexcept {
    return ceil_half(n - gain);
}

// ---------------------------------------------------------------------------
// ExactGainSolver
//
// Exact reformulation of the original objective:
//   cost(S) = |S| + ceil((n - |union(S)|) / 2)
// so for any target B,
//   cost(S) <= B  <=>  |union(S)| - 2|S| >= n - 2B.
//
// The search therefore maximizes exact gain = covered_points - 2 * chosen_lines,
// with no disjointness assumption.  It branches by binary include/exclude
// decisions on gain-positive lines instead of a point-partition tree.
// ---------------------------------------------------------------------------
class ExactGainSolver {
public:
    ExactGainSolver(const std::vector<HeavyLine>&            lines,
                    const std::vector<std::vector<int>>&     incidence)
        : lines_(lines)
        , incidence_(incidence)
        , worker_count_(std::max(1U, std::thread::hardware_concurrency())) {}

    [[nodiscard]] GainSolveResult solve_for_n(
        int n,
        const std::vector<int>& active_ids,
        const std::vector<int>& warm_start_line_ids,
        const std::vector<int>& greedy_seed_line_ids)
    {
        cur_n_  = n;
        words_  = words_for_n(n);
        best_gain_.store(0, std::memory_order_relaxed);
        best_choice_gain_ = 0;
        best_choice_.clear();

        ordered_point_lines_.assign(static_cast<std::size_t>(n), {});
        for (int p = 0; p < n; ++p) {
            auto& ordered = ordered_point_lines_[static_cast<std::size_t>(p)];
            ordered.reserve(incidence_[static_cast<std::size_t>(p)].size());
            for (int id : incidence_[static_cast<std::size_t>(p)]) {
                if (lines_[static_cast<std::size_t>(id)].activate_at <= n)
                    ordered.push_back(id);
            }
            std::stable_sort(ordered.begin(), ordered.end(), [&](int lhs, int rhs) {
                const auto sa = lines_[static_cast<std::size_t>(lhs)].points.size();
                const auto sb = lines_[static_cast<std::size_t>(rhs)].points.size();
                if (sa != sb) return sa > sb;
                return lhs < rhs;
            });
        }

        Task root = make_root_task(n, active_ids);
        seed_from_line_ids(root, warm_start_line_ids);
        seed_from_line_ids(root, greedy_seed_line_ids);

        const GreedyResult greedy = greedy_completion(root, true);
        if (greedy.gain > current_best_gain())
            submit_candidate(greedy.gain, greedy.lines);

        if (worker_count_ <= 1) {
            run_task(std::move(root));
        } else {
            std::vector<Task> frontier = build_frontier(std::move(root));
            if (frontier.empty()) {
                // Nothing left to search; incumbent came entirely from seeding.
            } else if (frontier.size() == 1) {
                run_task(std::move(frontier.front()));
            } else {
                solve_frontier_parallel(frontier);
            }
        }

        GainSolveResult out;
        out.best_gain = current_best_gain();
        out.chosen_line_ids = best_choice_;
        return out;
    }

private:
    struct Task {
        BitMask1024                   active;
        int                           active_count  = 0;
        int                           current_gain  = 0;
        std::vector<std::uint16_t>    line_cover;
        std::vector<unsigned char>    available;
        std::vector<int>              current_choice;
    };

    struct GreedyResult {
        int              gain = 0;
        std::vector<int> lines;
    };

    struct UndoLog {
        std::vector<int> disabled_lines;
        std::vector<int> decremented_lines;
        std::vector<int> newly_covered_points;
    };

    [[nodiscard]] Task make_root_task(
        int n, const std::vector<int>& active_ids) const
    {
        Task root;
        root.active = BitMask1024::prefix(n);
        root.active_count = n;
        root.line_cover.assign(lines_.size(), 0);
        root.available.assign(lines_.size(), 0);
        for (int id : active_ids) {
            const int cover = lines_[static_cast<std::size_t>(id)].mask.intersection_count(root.active, words_);
            root.line_cover[static_cast<std::size_t>(id)] = static_cast<std::uint16_t>(cover);
            root.available[static_cast<std::size_t>(id)] = 1;
        }
        return root;
    }

    [[nodiscard]] int current_best_gain() const noexcept {
        return best_gain_.load(std::memory_order_relaxed);
    }

    [[nodiscard]] int current_best_cost() const noexcept {
        return cost_from_gain(cur_n_, current_best_gain());
    }

    [[nodiscard]] int gain_needed_to_beat_best_cost() const noexcept {
        return cur_n_ - 2 * current_best_cost() + 2;
    }

    [[nodiscard]] bool gain_upper_bound_cannot_beat_best_cost(long double total_gain_ub) const noexcept {
        const int integer_gain_ub = static_cast<int>(std::floor(total_gain_ub + 1e-9L));
        return integer_gain_ub < gain_needed_to_beat_best_cost();
    }

    void submit_candidate(int gain, const std::vector<int>& choice) {
        int observed = current_best_gain();
        while (gain > observed &&
               !best_gain_.compare_exchange_weak(
                   observed, gain,
                   std::memory_order_relaxed,
                   std::memory_order_relaxed)) {}
        if (gain < current_best_gain()) return;
        std::lock_guard<std::mutex> lock(best_choice_mutex_);
        if (gain > best_choice_gain_) {
            best_choice_gain_ = gain;
            best_choice_ = choice;
        }
    }

    [[nodiscard]] int positive_gain(const Task& task, int id) const noexcept {
        if (!task.available[static_cast<std::size_t>(id)]) return 0;
        const int gain = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]) - 2;
        return gain > 0 ? gain : 0;
    }

    void disable_line(Task& task, UndoLog& undo, int id) const {
        unsigned char& available = task.available[static_cast<std::size_t>(id)];
        if (!available) return;
        available = 0;
        undo.disabled_lines.push_back(id);
    }

    void cover_point(Task& task, UndoLog& undo, int p) const {
        if (!task.active.test(p)) return;
        task.active.reset(p);
        --task.active_count;
        undo.newly_covered_points.push_back(p);
        for (int id : incidence_[static_cast<std::size_t>(p)]) {
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
            if (!task.available[static_cast<std::size_t>(id)]) continue;
            --task.line_cover[static_cast<std::size_t>(id)];
            undo.decremented_lines.push_back(id);
        }
    }

    void select_line(Task& task, UndoLog& undo, int id) const {
        disable_line(task, undo, id);
        for (int p : lines_[static_cast<std::size_t>(id)].points) {
            if (p >= cur_n_) break;
            cover_point(task, undo, p);
        }
    }

    void undo_to(Task& task, const UndoLog& undo) const {
        for (int id : undo.decremented_lines)
            ++task.line_cover[static_cast<std::size_t>(id)];
        for (int p : undo.newly_covered_points) {
            task.active.set(p);
            ++task.active_count;
        }
        for (int id : undo.disabled_lines)
            task.available[static_cast<std::size_t>(id)] = 1;
    }

    void select_line_without_undo(Task& task, int id) const {
        if (!task.available[static_cast<std::size_t>(id)]) return;
        task.available[static_cast<std::size_t>(id)] = 0;
        for (int p : lines_[static_cast<std::size_t>(id)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            task.active.reset(p);
            --task.active_count;
            for (int lid : incidence_[static_cast<std::size_t>(p)]) {
                if (lines_[static_cast<std::size_t>(lid)].activate_at > cur_n_) continue;
                if (!task.available[static_cast<std::size_t>(lid)]) continue;
                --task.line_cover[static_cast<std::size_t>(lid)];
            }
        }
    }

    void seed_from_line_ids(const Task& root, const std::vector<int>& line_ids) {
        if (line_ids.empty()) return;
        Task trial = root;
        int gain = 0;
        std::vector<int> chosen;
        chosen.reserve(line_ids.size());
        for (int id : line_ids) {
            if (id < 0 || id >= static_cast<int>(lines_.size())) continue;
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
            const int delta = positive_gain(trial, id);
            if (delta <= 0) continue;
            gain += delta;
            chosen.push_back(id);
            select_line_without_undo(trial, id);
        }
        if (gain > current_best_gain())
            submit_candidate(gain, chosen);
    }

    [[nodiscard]] int choose_branch_line(const Task& task) const
    {
        int best_line = -1;
        int best_gain = -1;
        long long best_pressure = -1;

        for (int id = 0; id < static_cast<int>(lines_.size()); ++id) {
            const int gain = positive_gain(task, id);
            if (gain <= 0) continue;

            long long pressure = 0;
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= cur_n_) break;
                if (!task.active.test(p)) continue;
                pressure += static_cast<long long>(
                    ordered_point_lines_[static_cast<std::size_t>(p)].size());
            }

            const bool better =
                (gain > best_gain) ||
                (gain == best_gain && pressure > best_pressure) ||
                (gain == best_gain && pressure == best_pressure && id < best_line);
            if (better) {
                best_line = id;
                best_gain = gain;
                best_pressure = pressure;
            }
        }
        return best_line;
    }

    [[nodiscard]] GreedyResult greedy_completion(const Task& task, bool keep_lines) const {
        using HeapItem = std::pair<int, int>;

        BitMask1024 temp_active = task.active;
        std::vector<unsigned char> temp_available = task.available;
        std::priority_queue<HeapItem> heap;
        for (int id = 0; id < static_cast<int>(lines_.size()); ++id) {
            if (!temp_available[static_cast<std::size_t>(id)]) continue;
            const int cover = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
            if (cover >= 3) heap.push({cover - 2, id});
        }

        GreedyResult result;
        while (!heap.empty()) {
            const int id = heap.top().second;
            heap.pop();
            if (!temp_available[static_cast<std::size_t>(id)]) continue;

            int current_gain = -2;
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= cur_n_) break;
                if (temp_active.test(p)) ++current_gain;
            }
            if (current_gain <= 0) {
                temp_available[static_cast<std::size_t>(id)] = 0;
                continue;
            }

            result.gain += current_gain;
            if (keep_lines) result.lines.push_back(id);
            temp_available[static_cast<std::size_t>(id)] = 0;
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= cur_n_) break;
                if (temp_active.test(p)) temp_active.reset(p);
            }
        }
        return result;
    }

    [[nodiscard]] long double lagrangian_upper_bound(const Task& task) const {
        std::vector<int> active_points;
        active_points.reserve(static_cast<std::size_t>(task.active_count));
        task.active.for_each_set_bit(words_, [&](int p) {
            bool seen = false;
            for (int id : incidence_[static_cast<std::size_t>(p)]) {
                if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
                if (positive_gain(task, id) > 0) {
                    seen = true;
                    break;
                }
            }
            if (seen) active_points.push_back(p);
        });
        if (active_points.empty())
            return static_cast<long double>(task.current_gain);

        std::vector<long double> y(static_cast<std::size_t>(cur_n_), 0.0L);
        for (int p : active_points) {
            long double best_density = 0.0L;
            for (int id : incidence_[static_cast<std::size_t>(p)]) {
                if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
                if (!task.available[static_cast<std::size_t>(id)]) continue;
                const int rem = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
                if (rem < 3) continue;
                const long double density =
                    static_cast<long double>(rem - 2) / static_cast<long double>(rem);
                best_density = std::max(best_density, density);
            }
            y[static_cast<std::size_t>(p)] = std::min<long double>(1.0L, best_density);
        }

        const int greedy_lb = greedy_completion(task, false).gain;
        auto eval = [&](const std::vector<long double>& yy,
                        std::vector<long double>* subgrad) -> long double {
            long double value = 0.0L;
            if (subgrad != nullptr)
                std::fill(subgrad->begin(), subgrad->end(), 0.0L);

            for (int p : active_points) {
                value += yy[static_cast<std::size_t>(p)];
                if (subgrad != nullptr)
                    (*subgrad)[static_cast<std::size_t>(p)] = 1.0L;
            }

            for (int id = 0; id < static_cast<int>(lines_.size()); ++id) {
                if (!task.available[static_cast<std::size_t>(id)]) continue;
                const int rem = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
                if (rem < 3) continue;

                long double sum_y = 0.0L;
                for (int p : lines_[static_cast<std::size_t>(id)].points) {
                    if (p >= cur_n_) break;
                    if (task.active.test(p))
                        sum_y += yy[static_cast<std::size_t>(p)];
                }
                const long double rc = static_cast<long double>(rem - 2) - sum_y;
                if (rc > 1e-18L) {
                    value += rc;
                    if (subgrad != nullptr) {
                        for (int p : lines_[static_cast<std::size_t>(id)].points) {
                            if (p >= cur_n_) break;
                            if (task.active.test(p))
                                (*subgrad)[static_cast<std::size_t>(p)] -= 1.0L;
                        }
                    }
                }
            }
            return value;
        };

        long double best = eval(y, nullptr);
        const long double initial_total = static_cast<long double>(task.current_gain) + best;
        if (gain_upper_bound_cannot_beat_best_cost(initial_total))
            return initial_total;

        std::vector<long double> subgrad(static_cast<std::size_t>(cur_n_), 0.0L);
        for (int iter = 0; iter < 36; ++iter) {
            const long double value = eval(y, &subgrad);
            best = std::min(best, value);

            long double norm2 = 0.0L;
            for (int p : active_points) {
                const long double g = subgrad[static_cast<std::size_t>(p)];
                norm2 += g * g;
            }
            if (norm2 <= 1e-30L) break;

            const long double incumbent_gap =
                std::max(0, gain_needed_to_beat_best_cost() - task.current_gain);
            const long double target = std::max<long double>(
                static_cast<long double>(greedy_lb), incumbent_gap);
            const long double gap = std::max<long double>(0.0L, value - target);
            if (gap <= 1e-15L) break;

            const long double step = (1.35L * gap) / norm2;
            for (int p : active_points) {
                long double next = y[static_cast<std::size_t>(p)]
                                 - step * subgrad[static_cast<std::size_t>(p)];
                if (next < 0.0L) next = 0.0L;
                else if (next > 1.0L) next = 1.0L;
                y[static_cast<std::size_t>(p)] = next;
            }
        }

        return static_cast<long double>(task.current_gain) + best;
    }

    void maybe_submit_greedy(const Task& task) {
        const GreedyResult greedy = greedy_completion(task, true);
        if (task.current_gain + greedy.gain <= current_best_gain()) return;

        std::vector<int> choice = task.current_choice;
        choice.insert(choice.end(), greedy.lines.begin(), greedy.lines.end());
        submit_candidate(task.current_gain + greedy.gain, choice);
    }

    [[nodiscard]] bool is_duplicate_active_coverage(
        const BitMask1024& rem, const std::vector<BitMask1024>& seen) const noexcept
    {
        for (const BitMask1024& prior : seen) {
            bool same = true;
            for (int wi = 0; wi < words_; ++wi) {
                if (prior.words[static_cast<std::size_t>(wi)] !=
                    rem.words[static_cast<std::size_t>(wi)]) {
                    same = false;
                    break;
                }
            }
            if (same) return true;
        }
        return false;
    }

    void dfs(Task& task) {
        const long double ub = lagrangian_upper_bound(task);
        if (gain_upper_bound_cannot_beat_best_cost(ub)) return;

        maybe_submit_greedy(task);

        const int branch_line = choose_branch_line(task);
        if (branch_line < 0) {
            submit_candidate(task.current_gain, task.current_choice);
            return;
        }

        {
            UndoLog undo;
            const int delta = positive_gain(task, branch_line);
            if (delta > 0) {
                select_line(task, undo, branch_line);
                task.current_gain += delta;
                task.current_choice.push_back(branch_line);
                dfs(task);
                task.current_choice.pop_back();
                task.current_gain -= delta;
            }
            undo_to(task, undo);
        }

        {
            UndoLog undo;
            disable_line(task, undo, branch_line);
            dfs(task);
            undo_to(task, undo);
        }
    }

    void run_task(Task task) {
        maybe_submit_greedy(task);
        dfs(task);
    }

    void expand_task(const Task& task, std::vector<Task>& children) const {
        children.clear();
        const int branch_line = choose_branch_line(task);
        if (branch_line < 0) return;

        children.reserve(2);

        Task include = task;
        const int delta = positive_gain(include, branch_line);
        if (delta > 0) {
            include.current_gain += delta;
            include.current_choice.push_back(branch_line);
            select_line_without_undo(include, branch_line);
            children.push_back(std::move(include));
        }

        Task exclude = task;
        exclude.available[static_cast<std::size_t>(branch_line)] = 0;
        children.push_back(std::move(exclude));
    }

    [[nodiscard]] std::vector<Task> build_frontier(Task root) {
        std::vector<Task> frontier;
        frontier.push_back(std::move(root));

        std::vector<Task> children;
        const std::size_t target =
            std::max<std::size_t>(1, static_cast<std::size_t>(worker_count_) * 2U);

        while (frontier.size() < target) {
            const auto it = std::max_element(
                frontier.begin(), frontier.end(),
                [](const Task& a, const Task& b) {
                    if (a.active_count != b.active_count) return a.active_count < b.active_count;
                    return a.current_gain > b.current_gain;
                });
            if (it == frontier.end()) break;

            Task task = std::move(*it);
            frontier.erase(it);

            const long double ub = lagrangian_upper_bound(task);
            if (gain_upper_bound_cannot_beat_best_cost(ub)) {
                if (frontier.empty()) break;
                continue;
            }

            maybe_submit_greedy(task);

            const int branch_line = choose_branch_line(task);
            if (branch_line < 0) {
                submit_candidate(task.current_gain, task.current_choice);
                if (frontier.empty()) break;
                continue;
            }

            expand_task(task, children);
            for (Task& child : children)
                frontier.push_back(std::move(child));
        }

        return frontier;
    }

    void solve_frontier_parallel(const std::vector<Task>& frontier) {
        const unsigned threads = std::min<unsigned>(
            worker_count_, static_cast<unsigned>(frontier.size()));
        if (threads <= 1) {
            run_task(frontier.front());
            return;
        }

        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(threads - 1));

        auto worker = [&]() {
            while (true) {
                const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= frontier.size()) break;
                run_task(frontier[idx]);
            }
        };

        for (unsigned i = 1; i < threads; ++i)
            workers.emplace_back(worker);
        worker();
        for (std::thread& thread : workers) thread.join();
    }

    const std::vector<HeavyLine>&           lines_;
    const std::vector<std::vector<int>>&    incidence_;
    unsigned                                worker_count_ = 1;
    int                                     cur_n_ = 0;
    int                                     words_ = 0;
    std::vector<std::vector<int>>           ordered_point_lines_;
    std::atomic<int>                        best_gain_{0};
    int                                     best_choice_gain_ = 0;
    std::vector<int>                        best_choice_;
    std::mutex                              best_choice_mutex_;
};

// ---------------------------------------------------------------------------
// MonolithicExactSolver
// Legacy direct-cost search family retained as an exact reference.
// The live solve_for_n() path currently uses ExactGainSolver above; this
// solver is kept for side-by-side experimentation, not for production runs.
// Incremental state: productive_count / size_freq are maintained via
// apply_delta, while per-depth state frames reuse storage across siblings.
// ---------------------------------------------------------------------------
class MonolithicExactSolver {
public:
    MonolithicExactSolver(const std::vector<HeavyLine>&            lines,
                          const std::vector<std::vector<int>>&     incidence)
        : lines_(lines)
        , incidence_(incidence)
        , worker_count_(std::max(1U, std::thread::hardware_concurrency())) {}

    [[nodiscard]] int solve_for_n(int n, const std::vector<int>& active_ids, int incumbent) {
        cur_n_  = n;
        words_  = words_for_n(n);
        best_.store(incumbent, std::memory_order_relaxed);
        StateKeyWithParityHash::used_words = words_;

        // Build the initial incremental state for the full active set.
        SearchTask root;
        root.active       = BitMask1024::prefix(n);
        root.active_count = n;
        root.dangling     = false;
        root.cost         = 0;
        reset_state(root.state);
        for (int id : active_ids) {
            const int c = lines_[id].mask.intersection_count(root.active, words_);
            root.state.line_cover[static_cast<std::size_t>(id)] = static_cast<std::uint16_t>(c);
            if (c >= 3) {
                root.state.size_freq[static_cast<std::size_t>(c)]++;
                for (int p : lines_[id].points) if (p < n) {
                    root.state.productive_count[static_cast<std::size_t>(p)]++;
                }
            }
        }

        if (worker_count_ <= 1) {
            WorkerContext ctx;
            ctx.state_stack.reserve(static_cast<std::size_t>(n + 1));
            run_task(root, ctx);
            return current_best();
        }

        std::vector<SearchTask> frontier = build_frontier(std::move(root));
        if (frontier.empty()) return current_best();
        if (frontier.size() == 1) {
            WorkerContext ctx;
            ctx.state_stack.reserve(static_cast<std::size_t>(n + 1));
            run_task(frontier.front(), ctx);
            return current_best();
        }
        return solve_frontier_parallel(frontier);
    }

private:
    struct SearchTask {
        BitMask1024      active;
        int              active_count = 0;
        bool             dangling     = false;
        int              cost         = 0;
        IncrementalState state;
    };

    struct WorkerContext {
        std::vector<IncrementalState> state_stack;
        std::unordered_map<StateKeyWithParity, int, StateKeyWithParityHash> memo;
    };

    struct DualBoundScratch {
        std::vector<int>           active_points;
        std::vector<int>           best_cover;
        std::vector<int>           productive_ids;
        std::vector<std::uint32_t> line_stamp;
        std::vector<long double>   slack;
        std::uint32_t              stamp = 0;

        void ensure(std::size_t line_count, int current_n) {
            if (best_cover.size() != static_cast<std::size_t>(current_n))
                best_cover.resize(static_cast<std::size_t>(current_n));
            if (line_stamp.size() != line_count) {
                line_stamp.assign(line_count, 0);
                slack.resize(line_count);
                stamp = 0;
            } else if (slack.size() != line_count) {
                slack.resize(line_count);
            }
        }

        [[nodiscard]] std::uint32_t next_stamp() {
            if (stamp == std::numeric_limits<std::uint32_t>::max()) {
                std::fill(line_stamp.begin(), line_stamp.end(), 0);
                stamp = 1;
                return stamp;
            }
            return ++stamp;
        }
    };

    IncrementalState& ensure_state_frame(WorkerContext& ctx, int depth) {
        if (depth >= static_cast<int>(ctx.state_stack.size()))
            ctx.state_stack.resize(static_cast<std::size_t>(depth + 1));
        IncrementalState& frame = ctx.state_stack[static_cast<std::size_t>(depth)];
        resize_state(frame);
        return frame;
    }

    void resize_state(IncrementalState& state) const {
        if (state.line_cover.size() != lines_.size())
            state.line_cover.resize(lines_.size());
        if (state.productive_count.size() != static_cast<std::size_t>(cur_n_))
            state.productive_count.resize(static_cast<std::size_t>(cur_n_));
        if (state.size_freq.size() != static_cast<std::size_t>(cur_n_ + 1))
            state.size_freq.resize(static_cast<std::size_t>(cur_n_ + 1));
    }

    void reset_state(IncrementalState& state) const {
        resize_state(state);
        std::fill(state.line_cover.begin(),       state.line_cover.end(),       0);
        std::fill(state.productive_count.begin(), state.productive_count.end(), 0);
        std::fill(state.size_freq.begin(),        state.size_freq.end(),        0);
    }

    static void copy_state(IncrementalState& dst, const IncrementalState& src) {
        std::copy(src.line_cover.begin(),       src.line_cover.end(),       dst.line_cover.begin());
        std::copy(src.productive_count.begin(), src.productive_count.end(), dst.productive_count.begin());
        std::copy(src.size_freq.begin(),        src.size_freq.end(),        dst.size_freq.begin());
    }

    [[nodiscard]] int current_best() const noexcept {
        return best_.load(std::memory_order_relaxed);
    }

    void tighten_best(int candidate) const noexcept {
        int observed = current_best();
        while (candidate < observed &&
               !best_.compare_exchange_weak(
                   observed, candidate,
                   std::memory_order_relaxed,
                   std::memory_order_relaxed)) {}
    }

    // Propagate the removal of all bits set in `removed` into state `s`.
    // For each removed point, decrement the cover of every line covering it.
    // When a line's cover drops below 3 it is no longer productive, so
    // decrement productive_count for its remaining points.
    void apply_delta(IncrementalState& s, const BitMask1024& removed) const {
        removed.for_each_set_bit(words_, [&](int p) {
            for (int id : incidence_[static_cast<std::size_t>(p)]) {
                if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
                std::uint16_t& lc = s.line_cover[static_cast<std::size_t>(id)];
                if (lc < 3) continue;
                s.size_freq[static_cast<std::size_t>(lc)]--;
                --lc;
                if (lc >= 3) {
                    s.size_freq[static_cast<std::size_t>(lc)]++;
                } else {
                    // Line is no longer productive; withdraw its contribution.
                    for (int lp : lines_[static_cast<std::size_t>(id)].points)
                        if (lp < cur_n_)
                            s.productive_count[static_cast<std::size_t>(lp)]--;
                }
            }
        });
    }

    static DualBoundScratch& dual_bound_scratch() {
        thread_local DualBoundScratch scratch;
        return scratch;
    }

    // Fractional lower bound based only on productive line cover sizes.
    // It is fast but ignores overlap structure between lines.
    [[nodiscard]] int size_histogram_lower_bound(
        const IncrementalState& state, int ac, bool d) const
    {
        const int pending = ac + static_cast<int>(d);
        int best = (pending + 1) / 2, covered = 0, chosen = 0, max_c = 2;
        for (int i = ac; i >= 3; --i)
            if (state.size_freq[static_cast<std::size_t>(i)] > 0) { max_c = i; break; }
        for (int i = ac; i >= 3; --i) {
            for (int j = 0; j < state.size_freq[static_cast<std::size_t>(i)]; ++j) {
                ++chosen;
                covered = std::min(ac, covered + i);
                best    = std::min(best, chosen + ceil_half(std::max(0, pending - covered)));
                if (covered == ac) return std::max(best, (pending + max_c - 1) / max_c);
            }
        }
        return std::max(best, (pending + max_c - 1) / max_c);
    }

    // Overlap-aware lower bound. It constructs dual-feasible point weights
    // for the LP relaxation
    //   min sum(line_vars) + 0.5 * sum(residual_point_vars)
    // whose dual is
    //   max sum(w_p)
    //   s.t. sum_{p in line} w_p <= 1 for every productive line,
    //        0 <= w_p <= 0.5 for every active point.
    // Any feasible dual value, plus 0.5 for an existing dangling point,
    // lower-bounds the exact remaining integer cost by weak duality.
    [[nodiscard]] int overlap_dual_lower_bound(
        const BitMask1024& active, const IncrementalState& state, bool dangling) const
    {
        DualBoundScratch& scratch = dual_bound_scratch();
        scratch.ensure(lines_.size(), cur_n_);

        scratch.active_points.clear();
        scratch.productive_ids.clear();
        const std::uint32_t stamp = scratch.next_stamp();

        active.for_each_set_bit(words_, [&](int p) {
            scratch.active_points.push_back(p);
            int best_cover = 0;
            for (int id : incidence_[static_cast<std::size_t>(p)]) {
                const int cover = state_line_cover(state, static_cast<std::size_t>(id));
                if (cover < 3) continue;
                best_cover = std::max(best_cover, cover);
                if (scratch.line_stamp[static_cast<std::size_t>(id)] != stamp) {
                    scratch.line_stamp[static_cast<std::size_t>(id)] = stamp;
                    scratch.productive_ids.push_back(id);
                }
            }
            scratch.best_cover[static_cast<std::size_t>(p)] = best_cover;
        });

        if (scratch.productive_ids.empty())
            return static_cast<int>(dangling);

        auto reset_slack = [&]() {
            for (int id : scratch.productive_ids)
                scratch.slack[static_cast<std::size_t>(id)] = 1.0L;
        };

        auto evaluate_order = [&](const auto& comp) -> long double {
            reset_slack();
            std::sort(scratch.active_points.begin(), scratch.active_points.end(), comp);

            long double total = 0.0L;
            for (int p : scratch.active_points) {
                long double weight = 0.5L;
                for (int id : incidence_[static_cast<std::size_t>(p)]) {
                    const int cover = state_line_cover(state, static_cast<std::size_t>(id));
                    if (cover < 3) continue;
                    weight = std::min(weight, scratch.slack[static_cast<std::size_t>(id)]);
                }
                total += weight;
                for (int id : incidence_[static_cast<std::size_t>(p)]) {
                    const int cover = state_line_cover(state, static_cast<std::size_t>(id));
                    if (cover < 3) continue;
                    scratch.slack[static_cast<std::size_t>(id)] -= weight;
                }
            }
            return total;
        };

        const auto by_best_cover_then_degree = [&](int a, int b) {
            const int ca = scratch.best_cover[static_cast<std::size_t>(a)];
            const int cb = scratch.best_cover[static_cast<std::size_t>(b)];
            if (ca != cb) return ca < cb;
            const int da = state.productive_count[static_cast<std::size_t>(a)];
            const int db = state.productive_count[static_cast<std::size_t>(b)];
            if (da != db) return da < db;
            return a < b;
        };
        const auto by_degree_then_best_cover = [&](int a, int b) {
            const int da = state.productive_count[static_cast<std::size_t>(a)];
            const int db = state.productive_count[static_cast<std::size_t>(b)];
            if (da != db) return da < db;
            const int ca = scratch.best_cover[static_cast<std::size_t>(a)];
            const int cb = scratch.best_cover[static_cast<std::size_t>(b)];
            if (ca != cb) return ca < cb;
            return a < b;
        };

        const long double dual_value = std::max(
            evaluate_order(by_best_cover_then_degree),
            evaluate_order(by_degree_then_best_cover));
        const long double remaining = dual_value + (dangling ? 0.5L : 0.0L);
        return static_cast<int>(std::ceil(remaining - 1e-12L));
    }

    [[nodiscard]] int lower_bound(
        const BitMask1024& active, const IncrementalState& state, int ac, bool dangling,
        int budget) const
    {
        const int coarse = size_histogram_lower_bound(state, ac, dangling);
        // If the cheap histogram bound already reaches the pruning threshold
        // (coarse >= budget means cost + coarse >= current_best), the dual
        // cannot change the outcome — skip it entirely.
        if (ac <= 0 || coarse >= budget) return coarse;
        return std::max(coarse, overlap_dual_lower_bound(active, state, dangling));
    }

    [[nodiscard]] bool normalize_task(
        SearchTask& task,
        std::unordered_map<StateKeyWithParity, int, StateKeyWithParityHash>& memo) const
    {
        tighten_best(task.cost + (task.active_count + static_cast<int>(task.dangling) + 1) / 2);
        if (task.cost >= current_best()) return false;

        while (task.active_count > 0) {
            BitMask1024 forced{};
            int  forced_count = 0;
            bool any = false;
            task.active.for_each_set_bit(words_, [&](int p) {
                if (task.state.productive_count[static_cast<std::size_t>(p)] > 0) any = true;
                else { forced.set(p); ++forced_count; }
            });
            if (!any) {
                tighten_best(task.cost + (task.active_count + static_cast<int>(task.dangling) + 1) / 2);
                return false;
            }
            if (forced_count == 0) {
                if (task.cost + lower_bound(
                        task.active, task.state, task.active_count, task.dangling,
                        current_best() - task.cost) >= current_best())
                    return false;
                break;
            }
            apply_delta(task.state, forced);
            task.active = task.active.and_not(forced);
            task.active_count -= forced_count;
            const int total = static_cast<int>(task.dangling) + forced_count;
            task.cost += total / 2;
            task.dangling = (total & 1) != 0;
            if (task.cost >= current_best()) return false;
        }
        if (task.active_count == 0) {
            tighten_best(task.cost + static_cast<int>(task.dangling));
            return false;
        }

        const StateKeyWithParity key{task.active.words, task.dangling};
        const auto [it, inserted] = memo.emplace(key, task.cost);
        if (!inserted) {
            if (it->second <= task.cost) return false;
            it->second = task.cost;
        }
        return true;
    }

    void collect_candidates(
        const BitMask1024& active, const IncrementalState& state,
        int& branch_point, std::vector<int>& candidates) const
    {
        branch_point = -1;
        int min_count = std::numeric_limits<int>::max();
        active.for_each_set_bit(words_, [&](int p) {
            const int cnt = state.productive_count[static_cast<std::size_t>(p)];
            if (cnt < min_count) {
                min_count    = cnt;
                branch_point = p;
            }
        });

        candidates.clear();
        if (branch_point < 0) return;
        for (int id : incidence_[static_cast<std::size_t>(branch_point)]) {
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
            if (state.line_cover[static_cast<std::size_t>(id)] >= 3)
                candidates.push_back(id);
        }
        std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
            const int ca = state.line_cover[static_cast<std::size_t>(a)];
            const int cb = state.line_cover[static_cast<std::size_t>(b)];
            if (ca != cb) return ca > cb;
            const auto sa = lines_[static_cast<std::size_t>(a)].points.size();
            const auto sb = lines_[static_cast<std::size_t>(b)].points.size();
            if (sa != sb) return sa > sb;
            return a < b;
        });
    }

    void expand_task(const SearchTask& task, std::vector<SearchTask>& children) const {
        int branch_point = -1;
        std::vector<int> candidates;
        collect_candidates(task.active, task.state, branch_point, candidates);
        children.clear();
        children.reserve(candidates.size() + 1);

        for (int id : candidates) {
            SearchTask child;
            child.active       = task.active.and_not(lines_[static_cast<std::size_t>(id)].mask);
            child.active_count = task.active_count -
                state_line_cover(task.state, static_cast<std::size_t>(id));
            child.dangling     = task.dangling;
            child.cost         = task.cost + 1;
            child.state        = task.state;
            const BitMask1024 rem = task.active & lines_[static_cast<std::size_t>(id)].mask;
            apply_delta(child.state, rem);
            children.push_back(std::move(child));
        }

        BitMask1024 skip_mask{};
        skip_mask.set(branch_point);
        SearchTask skip;
        skip.active       = task.active.and_not(skip_mask);
        skip.active_count = task.active_count - 1;
        skip.cost         = task.cost + (static_cast<int>(task.dangling) + 1) / 2;
        skip.dangling     = ((static_cast<int>(task.dangling) + 1) & 1) != 0;
        skip.state        = task.state;
        apply_delta(skip.state, skip_mask);
        children.push_back(std::move(skip));
    }

    [[nodiscard]] static int state_line_cover(const IncrementalState& state, std::size_t id) noexcept {
        return static_cast<int>(state.line_cover[id]);
    }

    [[nodiscard]] std::vector<SearchTask> build_frontier(SearchTask root) const {
        std::vector<SearchTask> frontier;
        frontier.push_back(std::move(root));

        std::unordered_map<StateKeyWithParity, int, StateKeyWithParityHash> split_memo;
        split_memo.reserve(1 << 12);

        std::vector<SearchTask> children;
        const std::size_t target =
            std::max<std::size_t>(1, static_cast<std::size_t>(worker_count_) * 2U);

        while (frontier.size() < target) {
            const auto it = std::max_element(
                frontier.begin(), frontier.end(),
                [](const SearchTask& a, const SearchTask& b) {
                    if (a.active_count != b.active_count) return a.active_count < b.active_count;
                    return a.cost > b.cost;
                });
            if (it == frontier.end()) break;

            SearchTask task = std::move(*it);
            frontier.erase(it);

            if (!normalize_task(task, split_memo)) {
                if (frontier.empty()) break;
                continue;
            }

            expand_task(task, children);
            for (SearchTask& child : children)
                frontier.push_back(std::move(child));
        }

        return frontier;
    }

    void run_task(const SearchTask& task, WorkerContext& ctx) {
        ctx.memo.clear();
        ctx.state_stack.clear();
        ctx.state_stack.reserve(static_cast<std::size_t>(cur_n_ + 1));
        IncrementalState& root = ensure_state_frame(ctx, 0);
        copy_state(root, task.state);
        dfs(task.active, task.active_count, task.dangling, task.cost, root, ctx, 0);
    }

    [[nodiscard]] int solve_frontier_parallel(const std::vector<SearchTask>& frontier) {
        const unsigned threads = std::min<unsigned>(
            worker_count_, static_cast<unsigned>(frontier.size()));
        if (threads <= 1) {
            WorkerContext ctx;
            ctx.state_stack.reserve(static_cast<std::size_t>(cur_n_ + 1));
            run_task(frontier.front(), ctx);
            return current_best();
        }

        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(threads - 1));

        auto worker = [&]() {
            WorkerContext ctx;
            ctx.state_stack.reserve(static_cast<std::size_t>(cur_n_ + 1));
            while (true) {
                const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= frontier.size()) break;
                run_task(frontier[idx], ctx);
            }
        };

        for (unsigned i = 1; i < threads; ++i)
            workers.emplace_back(worker);
        worker();
        for (std::thread& thread : workers) thread.join();
        return current_best();
    }

    void dfs(
        BitMask1024 active, int ac, bool d, int cost, IncrementalState& state,
        WorkerContext& ctx, int depth)
    {
        tighten_best(cost + (ac + static_cast<int>(d) + 1) / 2);
        if (cost >= current_best()) return;

        // Forced-propagation loop: points with no productive line are deferred
        // as isolated and paired cheaply.  Repeat until none remain or pruning
        // kills the branch.
        while (ac > 0) {
            BitMask1024 forced{};
            int  fc  = 0;
            bool any = false;
            active.for_each_set_bit(words_, [&](int p) {
                if (state.productive_count[static_cast<std::size_t>(p)] > 0) any = true;
                else { forced.set(p); ++fc; }
            });
            if (!any) {
                tighten_best(cost + (ac + static_cast<int>(d) + 1) / 2);
                return;
            }
            if (fc == 0) {
                if (cost + lower_bound(active, state, ac, d, current_best() - cost) >= current_best()) return;
                break;
            }
            apply_delta(state, forced);
            active = active.and_not(forced);
            ac -= fc;
            const int total = static_cast<int>(d) + fc;
            cost += total / 2;
            d     = (total & 1) != 0;
            if (cost >= current_best()) return;
        }
        if (ac == 0) { tighten_best(cost + static_cast<int>(d)); return; }

        // Memoize on (active mask, dangling parity): prune revisits at equal
        // or higher cost.
        const StateKeyWithParity key{active.words, d};
        const auto [it, inserted] = ctx.memo.emplace(key, cost);
        if (!inserted) {
            if (it->second <= cost) return;
            it->second = cost;
        }

        int bp = -1;
        std::vector<int> candidates;
        collect_candidates(active, state, bp, candidates);

        for (int id : candidates) {
            const int             cover = state.line_cover[static_cast<std::size_t>(id)];
            const BitMask1024     rem   = active & lines_[static_cast<std::size_t>(id)].mask;
            IncrementalState&     child = ensure_state_frame(ctx, depth + 1);
            copy_state(child, state);
            apply_delta(child, rem);
            dfs(active.and_not(lines_[static_cast<std::size_t>(id)].mask),
                ac - cover, d, cost + 1, child, ctx, depth + 1);
        }

        // Skip branch: defer bp as an isolated point.
        BitMask1024      skp{};
        skp.set(bp);
        IncrementalState& skip_s = ensure_state_frame(ctx, depth + 1);
        copy_state(skip_s, state);
        apply_delta(skip_s, skp);
        const int total = static_cast<int>(d) + 1;
        dfs(active.and_not(skp), ac - 1,
            (total & 1) != 0, cost + total / 2, skip_s, ctx, depth + 1);
    }

    const std::vector<HeavyLine>&        lines_;
    const std::vector<std::vector<int>>& incidence_;
    unsigned                worker_count_ = 1;
    int                     cur_n_ = 0;
    int                     words_ = 0;
    mutable std::atomic<int> best_{0};
};

// ---------------------------------------------------------------------------
// ComponentProfileSolver
// Used only when the heavy-line graph decomposes into truly independent
// sub-problems.  Memoises on active mask -> Summary2{even_cost, odd_cost}.
// Uses full-recompute scratch (acceptable for small sub-problems).
// ---------------------------------------------------------------------------
class ComponentProfileSolver {
public:
    ComponentProfileSolver(const std::vector<HeavyLine>& lines,
                           int                           current_n,
                           std::vector<int>              component_line_ids,
                           int                           line_cap)
        : lines_(lines)
        , current_n_(current_n)
        , used_words_(words_for_n(current_n))
        , line_ids_(std::move(component_line_ids))
        , line_cap_(line_cap)
        , inf_(line_cap + 1)
        , scratch_(static_cast<std::size_t>(current_n + 1))
    {
        for (int id : line_ids_)
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= current_n_) break;
                point_to_lines_[static_cast<std::size_t>(p)].push_back(id);
            }
    }

    [[nodiscard]] Summary2 solve(const BitMask1024& component_mask) {
        StateKeyHash::used_words = used_words_;
        memo_.clear();
        memo_.reserve(1 << 12);
        return solve_state(component_mask, 0);
    }

private:
    struct Branch { int line_id = -1, cover = 0; };

    // Per-depth scratch recomputed fully at each memoised node.
    struct FrameScratch {
        std::array<int, kBitCapacity>     productive_count{};
        std::array<int, kBitCapacity>     best_cover{};
        void clear(int n) {
            std::fill_n(productive_count.begin(), n, 0);
            std::fill_n(best_cover.begin(),       n, 0);
        }
    };

    const Summary2& solve_state(const BitMask1024& active, int depth) {
        const StateKey key{active.words};
        if (const auto it = memo_.find(key); it != memo_.end()) return it->second;

        const int ac = active.count(used_words_);
        Summary2  result{inf_, inf_};
        if (ac == 0) { result[0] = 0; return emplace(key, result); }

        FrameScratch& scratch = scratch_[static_cast<std::size_t>(depth)];
        scratch.clear(current_n_);

        bool        any_productive = false;
        BitMask1024 forced_defer{};
        int         forced_count = 0;

        for (int id : line_ids_) {
            const int c = lines_[static_cast<std::size_t>(id)]
                              .mask.intersection_count(active, used_words_);
            if (c < 3) continue;
            any_productive = true;
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= current_n_ || !active.test(p)) continue;
                ++scratch.productive_count[static_cast<std::size_t>(p)];
                scratch.best_cover[static_cast<std::size_t>(p)] =
                    std::max(scratch.best_cover[static_cast<std::size_t>(p)], c);
            }
        }

        if (!any_productive) {
            result[0] = ceil_half(ac);
            result[1] = ac / 2;
            return emplace(key, result);
        }

        active.for_each_set_bit(used_words_, [&](int p) {
            if (scratch.productive_count[static_cast<std::size_t>(p)] == 0) {
                forced_defer.set(p); ++forced_count;
            }
        });

        if (forced_count > 0) {
            const Summary2& child = solve_state(active.and_not(forced_defer), depth + 1);
            for (int par = 0; par < 2; ++par) {
                if (child[par] > line_cap_) continue;
                const int out_par = (forced_count + par) & 1;
                const int add     = (forced_count + par) / 2;
                result[out_par]   = std::min(result[out_par],
                                             capped_add(child[par], add, line_cap_));
            }
            return emplace(key, result);
        }

        // MRV with tiebreak on best_cover then point index.
        int bp  = -1;
        int mc  = std::numeric_limits<int>::max();
        int mbc = std::numeric_limits<int>::max();
        active.for_each_set_bit(used_words_, [&](int p) {
            const int cnt = scratch.productive_count[static_cast<std::size_t>(p)];
            const int bc  = scratch.best_cover[static_cast<std::size_t>(p)];
            if (cnt < mc || (cnt == mc && bc < mbc) ||
                (cnt == mc && bc == mbc && p > bp))
                { bp = p; mc = cnt; mbc = bc; }
        });

        std::vector<Branch> branches;
        branches.reserve(point_to_lines_[static_cast<std::size_t>(bp)].size());
        for (int id : point_to_lines_[static_cast<std::size_t>(bp)]) {
            const int c = lines_[static_cast<std::size_t>(id)]
                              .mask.intersection_count(active, used_words_);
            if (c >= 3) branches.push_back({id, c});
        }
        std::sort(branches.begin(), branches.end(), [&](const Branch& a, const Branch& b) {
            if (a.cover != b.cover) return a.cover > b.cover;
            const auto sa = lines_[static_cast<std::size_t>(a.line_id)].points.size();
            const auto sb = lines_[static_cast<std::size_t>(b.line_id)].points.size();
            if (sa != sb) return sa > sb;
            return a.line_id < b.line_id;
        });

        // Skip branch: defer bp as an isolated point.
        BitMask1024 skip_mask = active;
        skip_mask.reset(bp);
        const Summary2& skip = solve_state(skip_mask, depth + 1);
        for (int par = 0; par < 2; ++par) {
            if (skip[par] > line_cap_) continue;
            const int out_par = par ^ 1;
            const int add     = (par + 1) / 2;
            result[out_par]   = std::min(result[out_par],
                                         capped_add(skip[par], add, line_cap_));
        }

        for (const Branch& br : branches) {
            const BitMask1024& lmask = lines_[static_cast<std::size_t>(br.line_id)].mask;
            const Summary2& child = solve_state(active.and_not(lmask), depth + 1);
            for (int par = 0; par < 2; ++par)
                result[par] = std::min(result[par], capped_add(1, child[par], line_cap_));
        }

        return emplace(key, result);
    }

    const Summary2& emplace(const StateKey& key, const Summary2& val) {
        return memo_.emplace(key, val).first->second;
    }

    const std::vector<HeavyLine>& lines_;
    int              current_n_, used_words_;
    std::vector<int> line_ids_;
    int              line_cap_, inf_;
    std::array<std::vector<int>, kBitCapacity> point_to_lines_{};
    std::vector<FrameScratch> scratch_;
    std::unordered_map<StateKey, Summary2, StateKeyHash> memo_;
};

// ---------------------------------------------------------------------------
// Component decomposition
// ---------------------------------------------------------------------------
struct ComponentBuildResult {
    std::vector<ActiveComponent> components;
    int isolated_points = 0;
};

[[maybe_unused]] [[nodiscard]] ComponentBuildResult build_active_components(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
{
    ComponentBuildResult out;
    if (active_ids.empty()) { out.isolated_points = n; return out; }

    const int m = static_cast<int>(active_ids.size());
    std::vector<int> parent(static_cast<std::size_t>(m));
    std::iota(parent.begin(), parent.end(), 0);

    const auto find = [&](int x, auto&& self) -> int {
        if (parent[static_cast<std::size_t>(x)] == x) return x;
        parent[static_cast<std::size_t>(x)] = self(parent[static_cast<std::size_t>(x)], self);
        return parent[static_cast<std::size_t>(x)];
    };
    auto unite = [&](int a, int b) {
        a = find(a, find); b = find(b, find);
        if (a != b) parent[static_cast<std::size_t>(a)] = b;
    };

    std::vector<int> first_line(static_cast<std::size_t>(n), -1);
    for (int local = 0; local < m; ++local) {
        for (int p : lines[static_cast<std::size_t>(active_ids[static_cast<std::size_t>(local)])].points) {
            if (p >= n) break;
            int& owner = first_line[static_cast<std::size_t>(p)];
            if (owner == -1) owner = local;
            else unite(local, owner);
        }
    }

    std::unordered_map<int, int> root_to_comp;
    root_to_comp.reserve(static_cast<std::size_t>(m));
    for (int local = 0; local < m; ++local) {
        const int root         = find(local, find);
        const auto [it, ins]   = root_to_comp.emplace(root, static_cast<int>(out.components.size()));
        if (ins) out.components.push_back({});
        out.components[static_cast<std::size_t>(it->second)].line_ids
            .push_back(active_ids[static_cast<std::size_t>(local)]);
    }

    for (int p = 0; p < n; ++p) {
        const int owner = first_line[static_cast<std::size_t>(p)];
        if (owner == -1) { ++out.isolated_points; continue; }
        const int ci = root_to_comp[find(owner, find)];
        out.components[static_cast<std::size_t>(ci)].point_mask.set(p);
        ++out.components[static_cast<std::size_t>(ci)].point_count;
    }

    std::sort(out.components.begin(), out.components.end(),
              [](const ActiveComponent& a, const ActiveComponent& b) {
        if (a.point_count != b.point_count) return a.point_count > b.point_count;
        return a.line_ids.size() > b.line_ids.size();
    });
    return out;
}

[[nodiscard]] int solve_exact_with_components(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines,
    const std::vector<std::vector<int>>& incidence, int incumbent,
    const std::vector<int>& greedy_seed_line_ids, std::vector<int>& warm_start_solution)
{
    (void)incumbent;
    ExactGainSolver solver(lines, incidence);
    const GainSolveResult result =
        solver.solve_for_n(n, active_ids, warm_start_solution, greedy_seed_line_ids);
    warm_start_solution = result.chosen_line_ids;
    return (n - result.best_gain + 1) / 2;
}

} // namespace

// ---------------------------------------------------------------------------
// LiveDisplay — background thread that prints running progress to stderr.
// Only activates after a single N has been solving for more than 0.1 second,
// so it never clutters output for fast primes.  All writes go to stderr;
// stdout carries only the final per-N result lines.
//
// The main thread calls begin_n() just before each solve and clear() just
// after timing is captured.  Both are outside the timed region and consist
// of two cheap atomic stores — zero measurable impact on reported times.
// ---------------------------------------------------------------------------
class LiveDisplay {
public:
    using Clock = std::chrono::steady_clock;

    explicit LiveDisplay(Clock::time_point program_start)
        : program_start_(program_start)
    {
        thread_ = std::thread([this] { run(); });
    }

    ~LiveDisplay() {
        active_.store(false, std::memory_order_relaxed);
        thread_.join();
        clear_line();
    }

    // Called just before each solve starts (before t0 is captured).
    void begin_n(int n) noexcept {
        n_start_ns_.store(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - program_start_).count(),
            std::memory_order_relaxed);
        current_n_.store(n, std::memory_order_relaxed);
    }

    // Called after timing is captured, before printing the result line,
    // so the status line does not visually collide with the result.
    void clear() const noexcept { clear_line(); }

private:
    void clear_line() const {
        std::cerr << '\r' << std::string(78, ' ') << '\r' << std::flush;
    }

    void run() {
        while (active_.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
            if (!active_.load(std::memory_order_relaxed)) break;

            const int n = current_n_.load(std::memory_order_relaxed);
            if (n == 0) continue;

            const int64_t now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - program_start_).count();
            const double n_elapsed = static_cast<double>(
                now_ns - n_start_ns_.load(std::memory_order_relaxed)) * 1e-9;
            const double total = static_cast<double>(now_ns) * 1e-9;

            if (n_elapsed < 0.1) continue;  // silent during fast primes

            std::ostringstream oss;
            oss << "\r[running] N=" << n
                << "  N_time=" << std::fixed << std::setprecision(1) << n_elapsed << "s"
                << "  total="  << std::fixed << std::setprecision(1) << total    << "s  ";
            std::cerr << oss.str() << std::flush;
        }
    }

    Clock::time_point    program_start_;
    std::atomic<bool>    active_{true};
    std::atomic<int>     current_n_{0};
    std::atomic<int64_t> n_start_ns_{0};
    std::thread          thread_;
};

int main(int argc, char** argv) {
    int requested_n = kExecutionLimit;
    if (argc >= 2) {
        const int parsed = std::atoi(argv[1]);
        if (parsed > 0) requested_n = parsed;
    }
    requested_n = std::min({requested_n, kExecutionLimit, kBitCapacity});

    const std::vector<int>       primes = generate_primes(requested_n);
    const std::vector<HeavyLine> lines  = enumerate_heavy_lines(primes);

    std::vector<std::vector<int>> incidence( static_cast<std::size_t>(requested_n));
    std::vector<std::vector<int>> activation(static_cast<std::size_t>(requested_n + 1));
    for (int id = 0; id < static_cast<int>(lines.size()); ++id) {
        for (int p : lines[static_cast<std::size_t>(id)].points)
            if (p < requested_n)
                incidence[static_cast<std::size_t>(p)].push_back(id);
        if (lines[static_cast<std::size_t>(id)].activate_at <= requested_n)
            activation[static_cast<std::size_t>(lines[static_cast<std::size_t>(id)].activate_at)]
                .push_back(id);
    }

    std::vector<int> active_ids;
    active_ids.reserve(lines.size());

    const auto  program_start = std::chrono::steady_clock::now();
    LiveDisplay display(program_start);

    int  prev    = (kStartN > 1) ? kStartN - 1 : 0;
    bool stopped = false;
    std::vector<int> warm_start_solution;
    for (int n_init = 1; n_init < kStartN; ++n_init)
        for (int id : activation[static_cast<std::size_t>(n_init)])
            active_ids.push_back(id);
    for (int n = kStartN; n <= requested_n; ++n) {
        for (int id : activation[static_cast<std::size_t>(n)])
            active_ids.push_back(id);

        display.begin_n(n);                          // two atomic stores, before t0
        const auto t0 = std::chrono::steady_clock::now();

        const GreedyCoverSolution greedy = greedy_upper_bound_solution(n, active_ids, lines);
        int incumbent = std::min(n, prev + 1);
        incumbent     = std::min(incumbent, greedy.cost);

        const int answer = solve_exact_with_components(
            n, active_ids, lines, incidence, incumbent, greedy.line_ids, warm_start_solution);

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();

        display.clear();                             // wipe status line, after timing
        std::cout << "N="     << n
                  << " prime=" << primes[static_cast<std::size_t>(n - 1)]
                  << " lines=" << answer
                  << " time="  << format_seconds(elapsed)
                  << '\n';
        prev = answer;

        // Stop the sweep if this N exceeded the per-N time limit.
        // The check is after timing and output, so it never affects reported times.
        if (kPerNTimeLimitSeconds > 0.0 && elapsed > kPerNTimeLimitSeconds) {
            stopped = true;
            break;
        }
    }
    if (stopped)
        std::cout << "[stopped] N time limit of " << kPerNTimeLimitSeconds
                  << "s exceeded\n";
    return 0;
}