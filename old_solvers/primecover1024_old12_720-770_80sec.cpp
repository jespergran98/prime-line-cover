/*
 * primecover1024.cpp
 *
 * Exact minimum line cover solver for prime points (i, p_i).
 *
 * Live solver:
 *   - Enumerates every heavy line containing at least 3 points.
 *   - Solves the remaining exact-gain search with 1024-bit coverage masks,
 *     greedy seeding, a Lagrangian upper bound, and a shared atomic incumbent.
 *   - Uses a root cover-inequality lower bound before the branch-and-bound DFS.
 *   - Closes any residual uncovered points with exact 2-point / 1-point lines
 *     at cost ceil(residual / 2).
 *
 * Experiment Log:
 *   1. Adding an overlap-token ceiling inside lagrangian_upper_bound() was
 *      supposed to stop the DFS from giving the same contested uncovered
 *      points to too many future lines at once. In two full 770 runs the full
 *      version cut N=769 to 67880 nodes but slowed it to 27.707s / 27.198s,
 *      the near-prune gated version was slower still at 29.009s / 28.198s,
 *      and the lower-iteration follow-up hit the 30s per-N limit twice before
 *      it could finish N=769. The simplest explanation is that this per-node
 *      overlap accounting costs more than the extra pruning it buys on the
 *      hard tail, so future attempts should avoid broad contested-point
 *      exactification unless the same information can be reused cheaply or can
 *      be shown to lift the raw late-case lower bound itself.
 *   2. A branch-line overlap cap and two hotspot-only follow-ups were
 *      supposed to cut hard-tail dead branches by charging shared uncovered
 *      points more accurately inside the remaining-gain bound. In full 770
 *      runs they either slowed the sweep outright (85.572s and 83.836s) or
 *      only tied one run before losing the repeat (78.230s then 82.684s)
 *      against the 78.737s / 78.219s baseline, with N=769 staying pinned near
 *      21k nodes. The simplest explanation is that these extra point-overlap
 *      scans did not tighten the late raw bound enough to change the search
 *      tree, so future attempts should avoid adding more per-node overlap
 *      bookkeeping unless it can prune whole nodes much earlier or reuse the
 *      same work across many descendants.
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
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int    kBitCapacity          = 1024;
constexpr int    kBitWords             = kBitCapacity / 64;
constexpr int    kStartN               = 740;   // N to start from in the timed sweep; set to 0 to run from 1
constexpr int    kExecutionLimit       = 770;   // hard cap and default run-to N
// Maximum wall-clock seconds allowed per N. Set to 0.0 to disable (no limit). 60=1 min
constexpr double kPerNTimeLimitSeconds = 30;

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

// ---------------------------------------------------------------------------
// Core data structures
// ---------------------------------------------------------------------------
struct HeavyLine {
    BitMask1024      mask;
    std::vector<int> points;
    int              activate_at = 0;
};

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
    int              root_cost_lb = 0;
    int              raw_root_cost_lb = 0;
    long long        nodes = 0;
    std::vector<int> chosen_line_ids;
    bool             timed_out = false;
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
        , worker_count_(std::max(1U, std::thread::hardware_concurrency()))
        , total_incidence_(std::accumulate(
              lines.begin(), lines.end(), std::size_t{0},
              [](std::size_t total, const HeavyLine& line) {
                  return total + line.points.size();
              })) {}

    [[nodiscard]] GainSolveResult solve_for_n(
        int n,
        const std::vector<int>& active_ids,
        const std::vector<int>& warm_start_line_ids,
        const std::vector<int>& greedy_seed_line_ids,
        int root_cost_lb_floor)
    {
        cur_n_  = n;
        words_  = words_for_n(n);
        best_gain_.store(0, std::memory_order_relaxed);
        abort_.store(false, std::memory_order_relaxed);
        root_cost_lb_ = 0;
        raw_root_cost_lb_ = 0;
        root_cost_lb_floor_ = root_cost_lb_floor;
        root_lb_captured_ = false;
        node_count_.store(0, std::memory_order_relaxed);
        best_choice_gain_ = 0;
        best_choice_.clear();
        best_choice_.reserve(lines_.size());
        clear_seen_state_table();

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
        WorkerContext root_ctx;
        prepare_worker_context(root_ctx);
        seed_from_line_ids(root, warm_start_line_ids);
        seed_from_line_ids(root, greedy_seed_line_ids);

        const int greedy_gain = greedy_completion(root, root_ctx, true);
        if (greedy_gain > current_best_gain())
            submit_candidate(greedy_gain, root_ctx.greedy_lines);

        // Watchdog thread: sets abort_ once the per-N wall-clock limit is reached.
        // The flag is consumed inside lagrangian_upper_bound() — one relaxed load
        // amortised over all the expensive bound work — so there is zero overhead on
        // the DFS branch/select/undo hot path during normal operation.
        std::atomic<bool> watchdog_done{false};
        std::thread watchdog_thread;
        if (kPerNTimeLimitSeconds > 0.0) {
            watchdog_thread = std::thread([&]() {
                const auto deadline = std::chrono::steady_clock::now() +
                    std::chrono::duration<double>(kPerNTimeLimitSeconds);
                while (!watchdog_done.load(std::memory_order_relaxed)) {
                    const auto now = std::chrono::steady_clock::now();
                    if (now >= deadline) {
                        abort_.store(true, std::memory_order_relaxed);
                        return;
                    }
                    const auto remaining = deadline - now;
                    const auto slice = std::chrono::duration_cast<
                        std::chrono::steady_clock::duration>(
                            std::chrono::milliseconds(50));
                    std::this_thread::sleep_for(
                        remaining < slice ? remaining : slice);
                }
            });
        }

        if (worker_count_ <= 1) {
            run_task(std::move(root), root_ctx);
        } else {
            std::vector<Task> frontier = build_frontier(std::move(root), root_ctx);
            if (frontier.empty()) {
                // Nothing left to search; incumbent came entirely from seeding.
            } else if (frontier.size() == 1) {
                run_task(std::move(frontier.front()), root_ctx);
            } else {
                solve_frontier_parallel(frontier);
            }
        }

        // Signal the watchdog to exit and wait for it before we touch any locals.
        watchdog_done.store(true, std::memory_order_relaxed);
        if (watchdog_thread.joinable()) watchdog_thread.join();

        GainSolveResult out;
        out.best_gain     = current_best_gain();
        out.root_cost_lb = root_cost_lb_;
        out.raw_root_cost_lb = raw_root_cost_lb_;
        out.nodes        = node_count_.load(std::memory_order_relaxed);
        out.chosen_line_ids = best_choice_;
        out.timed_out     = abort_.load(std::memory_order_relaxed);
        return out;
    }

private:
    using HeapItem = std::pair<int, int>;

    struct Task {
        BitMask1024                   active;
        int                           active_count  = 0;
        int                           current_gain  = 0;
        std::vector<std::uint16_t>    line_cover;
        std::vector<unsigned char>    available;
        std::vector<int>              current_choice;
        // Unavailable line IDs along this path.  Duplicate-state pruning only
        // keeps the subset that is still productive under the current active
        // mask, because once a blocked line drops below gain-positive coverage
        // it can never matter again.
        std::vector<int>              blocked_line_ids;
        // Compact set of productive line IDs: available[id]==1 && line_cover[id]>=3.
        // productive_pos[id] = index in productive_ids, or -1 if not in set.
        std::vector<int>              productive_ids;
        std::vector<int>              productive_pos;
    };

    struct UndoMark {
        std::size_t disabled_lines = 0;
        std::size_t decremented_lines = 0;
        std::size_t newly_covered_points = 0;
        std::size_t removed_from_productive = 0;
        std::size_t blocked_lines = 0;
    };

    struct StateKey {
        BitMask1024      active;
        std::vector<int> blocked_productive_ids;

        bool operator==(const StateKey& other) const noexcept {
            return active.words == other.active.words &&
                   blocked_productive_ids == other.blocked_productive_ids;
        }
    };

    struct StateKeyHash {
        static std::uint64_t mix(std::uint64_t x) noexcept {
            x += 0x9e3779b97f4a7c15ULL;
            x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        }

        std::size_t operator()(const StateKey& key) const noexcept {
            std::uint64_t h = 0x9e3779b97f4a7c15ULL;
            for (std::uint64_t word : key.active.words)
                h ^= mix(word + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            for (int id : key.blocked_productive_ids)
                h ^= mix(static_cast<std::uint64_t>(id) +
                         0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            return static_cast<std::size_t>(h);
        }
    };

    using StateBestMap = std::unordered_map<StateKey, int, StateKeyHash>;

    struct SeenStateShard {
        std::mutex   mutex;
        StateBestMap best_gain;
    };

    struct LagrangianScratch {
        std::vector<int>         active_points;
        std::vector<long double> y;
        std::vector<long double> best_y;
        std::vector<long double> subgrad;
        std::vector<long double> line_slack;
        bool                     dual_ready = false;

        void ensure(std::size_t line_count, int current_n) {
            const std::size_t n = static_cast<std::size_t>(current_n);
            if (active_points.capacity() < n)
                active_points.reserve(n);
            if (y.size() != n)
                y.resize(n);
            if (best_y.size() != n)
                best_y.resize(n);
            if (subgrad.size() != n)
                subgrad.resize(n);
            if (line_slack.size() != line_count)
                line_slack.resize(line_count);
        }
    };

    struct WorkerContext {
        LagrangianScratch          lagrangian;
        std::vector<HeapItem>      greedy_heap;
        std::vector<unsigned char> greedy_available;
        std::vector<int>           greedy_lines;
        std::vector<int>           candidate_choice;
        std::vector<int>           disabled_lines;
        std::vector<int>           decremented_lines;
        std::vector<int>           newly_covered_points;
        std::vector<int>           removed_from_productive;
        std::vector<int>           state_blocked_ids;
        std::vector<int>           saved_dual_points;
        std::vector<long double>   saved_dual_values;
    };

    [[nodiscard]] Task make_root_task(
        int n, const std::vector<int>& active_ids) const
    {
        Task root;
        root.active = BitMask1024::prefix(n);
        root.active_count = n;
        root.line_cover.assign(lines_.size(), 0);
        root.available.assign(lines_.size(), 0);
        root.productive_pos.assign(lines_.size(), -1);
        root.current_choice.reserve(lines_.size());
        root.productive_ids.reserve(active_ids.size());
        for (int id : active_ids) {
            const int cover = lines_[static_cast<std::size_t>(id)].mask.intersection_count(root.active, words_);
            root.line_cover[static_cast<std::size_t>(id)] = static_cast<std::uint16_t>(cover);
            root.available[static_cast<std::size_t>(id)] = 1;
            if (cover >= 3) {
                root.productive_pos[static_cast<std::size_t>(id)] =
                    static_cast<int>(root.productive_ids.size());
                root.productive_ids.push_back(id);
            }
        }
        return root;
    }

    void ensure_task_search_buffers(Task& task) const {
        if (task.current_choice.capacity() < lines_.size())
            task.current_choice.reserve(lines_.size());
        if (task.blocked_line_ids.capacity() < 64)
            task.blocked_line_ids.reserve(64);
    }

    void prepare_worker_context(WorkerContext& ctx) const {
        ctx.lagrangian.ensure(lines_.size(), cur_n_);

        if (ctx.greedy_available.size() != lines_.size())
            ctx.greedy_available.resize(lines_.size());
        if (ctx.greedy_heap.capacity() < lines_.size())
            ctx.greedy_heap.reserve(lines_.size());
        ctx.greedy_heap.clear();

        if (ctx.greedy_lines.capacity() < lines_.size())
            ctx.greedy_lines.reserve(lines_.size());
        ctx.greedy_lines.clear();

        if (ctx.candidate_choice.capacity() < lines_.size())
            ctx.candidate_choice.reserve(lines_.size());
        ctx.candidate_choice.clear();

        if (ctx.disabled_lines.capacity() < lines_.size())
            ctx.disabled_lines.reserve(lines_.size());
        ctx.disabled_lines.clear();

        if (ctx.decremented_lines.capacity() < total_incidence_)
            ctx.decremented_lines.reserve(total_incidence_);
        ctx.decremented_lines.clear();

        if (ctx.newly_covered_points.capacity() < static_cast<std::size_t>(cur_n_))
            ctx.newly_covered_points.reserve(static_cast<std::size_t>(cur_n_));
        ctx.newly_covered_points.clear();

        if (ctx.removed_from_productive.capacity() < lines_.size())
            ctx.removed_from_productive.reserve(lines_.size());
        ctx.removed_from_productive.clear();

        if (ctx.state_blocked_ids.capacity() < 64)
            ctx.state_blocked_ids.reserve(64);
        ctx.state_blocked_ids.clear();

        if (ctx.saved_dual_points.capacity() < static_cast<std::size_t>(cur_n_))
            ctx.saved_dual_points.reserve(static_cast<std::size_t>(cur_n_));
        ctx.saved_dual_points.clear();

        if (ctx.saved_dual_values.capacity() < static_cast<std::size_t>(cur_n_))
            ctx.saved_dual_values.reserve(static_cast<std::size_t>(cur_n_));
        ctx.saved_dual_values.clear();

        ctx.lagrangian.dual_ready = false;
    }

    [[nodiscard]] int current_best_gain() const noexcept {
        return best_gain_.load(std::memory_order_relaxed);
    }

    [[nodiscard]] int current_best_cost() const noexcept {
        return cost_from_gain(cur_n_, current_best_gain());
    }

    [[nodiscard]] bool root_interval_closed() const noexcept {
        return root_lb_captured_ && current_best_cost() <= root_cost_lb_;
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

    [[nodiscard]] static UndoMark mark_undo(const Task& task, const WorkerContext& ctx) noexcept {
        return {
            ctx.disabled_lines.size(),
            ctx.decremented_lines.size(),
            ctx.newly_covered_points.size(),
            ctx.removed_from_productive.size(),
            task.blocked_line_ids.size()
        };
    }

    [[nodiscard]] std::size_t save_current_dual(WorkerContext& ctx) const {
        const std::size_t mark = ctx.saved_dual_points.size();
        if (!ctx.lagrangian.dual_ready) return mark;

        const auto& active_points = ctx.lagrangian.active_points;
        const auto& best_y = ctx.lagrangian.best_y;
        for (int p : active_points) {
            ctx.saved_dual_points.push_back(p);
            ctx.saved_dual_values.push_back(best_y[static_cast<std::size_t>(p)]);
        }
        return mark;
    }

    void restore_saved_dual(WorkerContext& ctx, std::size_t mark) const {
        for (std::size_t i = mark; i < ctx.saved_dual_points.size(); ++i) {
            const int p = ctx.saved_dual_points[i];
            ctx.lagrangian.best_y[static_cast<std::size_t>(p)] =
                ctx.saved_dual_values[i];
        }
        ctx.saved_dual_points.resize(mark);
        ctx.saved_dual_values.resize(mark);
        ctx.lagrangian.dual_ready = true;
    }

    // O(1) swap-to-back removal from productive_ids; records in the worker stack for reversal.
    void remove_productive(Task& task, WorkerContext& ctx, int id) const {
        int& pos = task.productive_pos[static_cast<std::size_t>(id)];
        if (pos < 0) return;
        const int last_id = task.productive_ids.back();
        task.productive_ids[static_cast<std::size_t>(pos)] = last_id;
        task.productive_pos[static_cast<std::size_t>(last_id)] = pos;
        task.productive_ids.pop_back();
        pos = -1;
        ctx.removed_from_productive.push_back(id);
    }

    // Called during undo to re-insert a line that was removed from the productive set.
    void reinsert_productive(Task& task, int id) const {
        if (task.productive_pos[static_cast<std::size_t>(id)] >= 0) return;
        task.productive_pos[static_cast<std::size_t>(id)] =
            static_cast<int>(task.productive_ids.size());
        task.productive_ids.push_back(id);
    }

    void disable_line_without_undo(Task& task, int id) const {
        unsigned char& available = task.available[static_cast<std::size_t>(id)];
        if (!available) return;
        available = 0;
        task.blocked_line_ids.push_back(id);
        int& pos = task.productive_pos[static_cast<std::size_t>(id)];
        if (pos >= 0) {
            const int last_id = task.productive_ids.back();
            task.productive_ids[static_cast<std::size_t>(pos)] = last_id;
            task.productive_pos[static_cast<std::size_t>(last_id)] = pos;
            task.productive_ids.pop_back();
            pos = -1;
        }
    }

    void disable_line(Task& task, WorkerContext& ctx, int id) const {
        unsigned char& available = task.available[static_cast<std::size_t>(id)];
        if (!available) return;
        available = 0;
        task.blocked_line_ids.push_back(id);
        ctx.disabled_lines.push_back(id);
        remove_productive(task, ctx, id);  // unavailable => non-productive
    }

    void cover_point(Task& task, WorkerContext& ctx, int p) const {
        if (!task.active.test(p)) return;
        task.active.reset(p);
        --task.active_count;
        ctx.newly_covered_points.push_back(p);
        for (int id : incidence_[static_cast<std::size_t>(p)]) {
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
            if (!task.available[static_cast<std::size_t>(id)]) continue;
            const uint16_t old_cover = task.line_cover[static_cast<std::size_t>(id)];
            --task.line_cover[static_cast<std::size_t>(id)];
            ctx.decremented_lines.push_back(id);
            if (old_cover == 3)   // cover 3->2: no longer productive
                remove_productive(task, ctx, id);
        }
    }

    void select_line(Task& task, WorkerContext& ctx, int id) const {
        disable_line(task, ctx, id);
        for (int p : lines_[static_cast<std::size_t>(id)].points) {
            if (p >= cur_n_) break;
            cover_point(task, ctx, p);
        }
    }

    void undo_to(Task& task, WorkerContext& ctx, UndoMark mark) const {
        // Restore productive set first (reverse of removal order), then restore
        // line_cover and available so the invariant is consistent throughout.
        for (std::size_t i = mark.removed_from_productive;
             i < ctx.removed_from_productive.size(); ++i)
            reinsert_productive(task, ctx.removed_from_productive[i]);
        for (std::size_t i = mark.decremented_lines;
             i < ctx.decremented_lines.size(); ++i)
            ++task.line_cover[static_cast<std::size_t>(ctx.decremented_lines[i])];
        for (std::size_t i = mark.newly_covered_points;
             i < ctx.newly_covered_points.size(); ++i) {
            const int p = ctx.newly_covered_points[i];
            task.active.set(p);
            ++task.active_count;
        }
        for (std::size_t i = mark.disabled_lines; i < ctx.disabled_lines.size(); ++i)
            task.available[static_cast<std::size_t>(ctx.disabled_lines[i])] = 1;
        task.blocked_line_ids.resize(mark.blocked_lines);

        ctx.removed_from_productive.resize(mark.removed_from_productive);
        ctx.decremented_lines.resize(mark.decremented_lines);
        ctx.newly_covered_points.resize(mark.newly_covered_points);
        ctx.disabled_lines.resize(mark.disabled_lines);
    }

    void select_line_without_undo(Task& task, int id) const {
        if (!task.available[static_cast<std::size_t>(id)]) return;
        disable_line_without_undo(task, id);

        for (int p : lines_[static_cast<std::size_t>(id)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            task.active.reset(p);
            --task.active_count;
            for (int lid : incidence_[static_cast<std::size_t>(p)]) {
                if (lines_[static_cast<std::size_t>(lid)].activate_at > cur_n_) continue;
                if (!task.available[static_cast<std::size_t>(lid)]) continue;
                const uint16_t old_cover = task.line_cover[static_cast<std::size_t>(lid)];
                --task.line_cover[static_cast<std::size_t>(lid)];
                if (old_cover == 3) {   // cover 3->2: drop from productive set
                    int& pos2 = task.productive_pos[static_cast<std::size_t>(lid)];
                    if (pos2 >= 0) {
                        const int last2 = task.productive_ids.back();
                        task.productive_ids[static_cast<std::size_t>(pos2)] = last2;
                        task.productive_pos[static_cast<std::size_t>(last2)] = pos2;
                        task.productive_ids.pop_back();
                        pos2 = -1;
                    }
                }
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

        // Iterate only over productive lines: available && cover >= 3 guaranteed by invariant.
        for (int id : task.productive_ids) {
            const int gain = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]) - 2;
            // gain > 0 guaranteed since cover >= 3 implies cover - 2 >= 1

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

    [[nodiscard]] int greedy_completion(
        const Task& task, WorkerContext& ctx, bool keep_lines) const
    {
        BitMask1024 temp_active = task.active;
        std::copy(task.available.begin(), task.available.end(), ctx.greedy_available.begin());

        auto& heap = ctx.greedy_heap;
        heap.clear();
        for (int id : task.productive_ids) {
            const int cover = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
            heap.push_back({cover - 2, id});
        }
        std::make_heap(heap.begin(), heap.end());

        auto& chosen = ctx.greedy_lines;
        chosen.clear();

        int gain = 0;
        while (!heap.empty()) {
            std::pop_heap(heap.begin(), heap.end());
            const int id = heap.back().second;
            heap.pop_back();
            if (!ctx.greedy_available[static_cast<std::size_t>(id)]) continue;

            int current_gain = -2;
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= cur_n_) break;
                if (temp_active.test(p)) ++current_gain;
            }
            if (current_gain <= 0) {
                ctx.greedy_available[static_cast<std::size_t>(id)] = 0;
                continue;
            }

            gain += current_gain;
            if (keep_lines) chosen.push_back(id);
            ctx.greedy_available[static_cast<std::size_t>(id)] = 0;
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= cur_n_) break;
                if (temp_active.test(p)) temp_active.reset(p);
            }
        }
        return gain;
    }

    [[nodiscard]] long double lagrangian_upper_bound(
        const Task& task, WorkerContext& ctx, int max_iters) const
    {
        // Zero-overhead cancellation: one relaxed atomic load, amortised over all the
        // expensive work below (greedy_completion + 36 subgradient iters).  Returning
        // -1 makes gain_upper_bound_cannot_beat_best_cost() prune this node immediately
        // in dfs() without touching any other hot-path code.
        if (abort_.load(std::memory_order_relaxed))
            return -1.0L;

        auto& scratch = ctx.lagrangian;
        auto& active_points = scratch.active_points;
        auto& y = scratch.y;
        auto& best_y = scratch.best_y;
        auto& subgrad = scratch.subgrad;
        auto& line_slack = scratch.line_slack;

        active_points.clear();
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
        if (active_points.empty()) {
            scratch.dual_ready = false;
            const long double total = static_cast<long double>(task.current_gain);
            return total;
        }

        if (scratch.dual_ready) {
            for (int p : active_points) {
                const long double warm = best_y[static_cast<std::size_t>(p)];
                y[static_cast<std::size_t>(p)] = std::clamp(warm, 0.0L, 1.0L);
            }
        } else {
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
        }

        auto eval = [&](const std::vector<long double>& yy,
                        std::vector<long double>* subgrad) -> long double {
            long double value = 0.0L;

            for (int p : active_points) {
                value += yy[static_cast<std::size_t>(p)];
                if (subgrad != nullptr)
                    (*subgrad)[static_cast<std::size_t>(p)] = 1.0L;
            }

            for (int id : task.productive_ids) {
                const int rem = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
                // available && rem >= 3 guaranteed by productive invariant — no guards needed.
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
        for (int p : active_points)
            best_y[static_cast<std::size_t>(p)] = y[static_cast<std::size_t>(p)];
        const long double initial_total = static_cast<long double>(task.current_gain) + best;
        scratch.dual_ready = true;
        if (gain_upper_bound_cannot_beat_best_cost(initial_total))
            return initial_total;

        const int greedy_lb = greedy_completion(task, ctx, false);
        for (int iter = 0; iter < max_iters; ++iter) {
            const long double value = eval(y, &subgrad);
            if (value + 1e-18L < best) {
                best = value;
                for (int p : active_points)
                    best_y[static_cast<std::size_t>(p)] = y[static_cast<std::size_t>(p)];
            }

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

        constexpr int kCoordinatePolishSweeps = 2;
        if (kCoordinatePolishSweeps > 0) {
            for (int id : task.productive_ids) {
                const int rem = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
                long double slack = static_cast<long double>(rem - 2);
                for (int p : lines_[static_cast<std::size_t>(id)].points) {
                    if (p >= cur_n_) break;
                    if (task.active.test(p))
                        slack -= best_y[static_cast<std::size_t>(p)];
                }
                line_slack[static_cast<std::size_t>(id)] = slack;
            }

            auto slack_objective = [&]() -> long double {
                long double value = 0.0L;
                for (int p : active_points)
                    value += best_y[static_cast<std::size_t>(p)];
                for (int id : task.productive_ids)
                    if (line_slack[static_cast<std::size_t>(id)] > 1e-18L)
                        value += line_slack[static_cast<std::size_t>(id)];
                return value;
            };

            long double polished_best = best;
            for (int sweep = 0; sweep < kCoordinatePolishSweeps; ++sweep) {
                bool changed = false;
                for (int p : active_points) {
                    const long double old = best_y[static_cast<std::size_t>(p)];
                    long double largest = 0.0L;
                    long double second  = 0.0L;
                    for (int id : incidence_[static_cast<std::size_t>(p)]) {
                        if (!task.available[static_cast<std::size_t>(id)]) continue;
                        const int rem = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
                        if (rem < 3) continue;

                        const long double adjusted =
                            line_slack[static_cast<std::size_t>(id)] + old;
                        if (adjusted <= 1e-18L) continue;
                        if (adjusted >= largest) {
                            second = largest;
                            largest = adjusted;
                        } else if (adjusted > second) {
                            second = adjusted;
                        }
                    }

                    const long double next = std::min<long double>(1.0L, second);
                    if (std::abs(next - old) <= 1e-15L) continue;
                    changed = true;
                    best_y[static_cast<std::size_t>(p)] = next;
                    const long double delta = next - old;
                    for (int id : incidence_[static_cast<std::size_t>(p)]) {
                        if (!task.available[static_cast<std::size_t>(id)]) continue;
                        const int rem = static_cast<int>(task.line_cover[static_cast<std::size_t>(id)]);
                        if (rem < 3) continue;
                        line_slack[static_cast<std::size_t>(id)] -= delta;
                    }
                }

                polished_best = std::min(polished_best, slack_objective());
                if (!changed) break;
            }
            best = polished_best;
        }

        return static_cast<long double>(task.current_gain) + best;
    }

    [[nodiscard]] long double lagrangian_upper_bound(
        const Task& task, WorkerContext& ctx) const
    {
        return lagrangian_upper_bound(task, ctx, 64); //optimal value ranges found around 62-66 and 112-116. 64 was the best single value in testing.
    }

    void capture_root_cost_lb(long double gain_ub) {
        if (root_lb_captured_) return;
        raw_root_cost_lb_ = cost_from_gain(
            cur_n_, static_cast<int>(std::floor(gain_ub + 1e-9L)));
        root_cost_lb_ = std::max(raw_root_cost_lb_, root_cost_lb_floor_);
        root_lb_captured_ = true;
    }

    void maybe_submit_greedy(const Task& task, WorkerContext& ctx) {
        const int greedy_gain = greedy_completion(task, ctx, true);
        if (task.current_gain + greedy_gain <= current_best_gain()) return;

        ctx.candidate_choice = task.current_choice;
        ctx.candidate_choice.insert(
            ctx.candidate_choice.end(),
            ctx.greedy_lines.begin(), ctx.greedy_lines.end());
        submit_candidate(task.current_gain + greedy_gain, ctx.candidate_choice);
    }

    [[nodiscard]] long double exclude_branch_upper_bound(
        const Task& task, const WorkerContext& ctx,
        int branch_line, long double node_gain_ub) const
    {
        const auto& best_y = ctx.lagrangian.best_y;
        const auto& line_slack = ctx.lagrangian.line_slack;

        long double capped = node_gain_ub;
        capped -= std::max(0.0L, line_slack[static_cast<std::size_t>(branch_line)]);

        for (int p : lines_[static_cast<std::size_t>(branch_line)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;

            const long double residual =
                1.0L - best_y[static_cast<std::size_t>(p)];
            if (residual <= 1e-18L) continue;

            long double requiring_sum = 0.0L;
            long double requiring_best = 0.0L;
            for (int id : incidence_[static_cast<std::size_t>(p)]) {
                if (id == branch_line) continue;
                if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
                if (!task.available[static_cast<std::size_t>(id)]) continue;
                if (task.line_cover[static_cast<std::size_t>(id)] < 3) continue;

                const long double rc = line_slack[static_cast<std::size_t>(id)];
                if (rc <= 1e-18L) continue;
                if (rc > residual + 1e-18L) continue;

                requiring_sum += rc;
                requiring_best = std::max(requiring_best, rc);
            }
            capped -= (requiring_sum - requiring_best);
        }

        return std::max<long double>(
            static_cast<long double>(task.current_gain), capped);
    }

    void clear_seen_state_table() {
        for (SeenStateShard& shard : seen_state_shards_)
            shard.best_gain.clear();
    }

    void collect_blocked_productive_lines(
        const Task& task, std::vector<int>& out) const
    {
        out.clear();
        if (task.blocked_line_ids.empty()) return;

        for (int id : task.blocked_line_ids) {
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
            const int cover =
                lines_[static_cast<std::size_t>(id)].mask.intersection_count(task.active, words_);
            if (cover >= 3)
                out.push_back(id);
        }
        if (out.size() > 1)
            std::sort(out.begin(), out.end());
    }

    [[nodiscard]] StateKey make_state_key(const Task& task, WorkerContext& ctx) const {
        collect_blocked_productive_lines(task, ctx.state_blocked_ids);
        StateKey key;
        key.active = task.active;
        key.blocked_productive_ids = ctx.state_blocked_ids;
        return key;
    }

    [[nodiscard]] bool seen_search_state(const Task& task, WorkerContext& ctx) {
        StateKey key = make_state_key(task, ctx);
        const std::size_t hash = StateKeyHash{}(key);
        SeenStateShard& shard =
            seen_state_shards_[hash % kSeenStateShardCount];
        std::lock_guard<std::mutex> lock(shard.mutex);
        const auto it = shard.best_gain.find(key);
        if (it != shard.best_gain.end()) {
            if (it->second >= task.current_gain) return true;
            it->second = task.current_gain;
            return false;
        }
        shard.best_gain.emplace(std::move(key), task.current_gain);
        return false;
    }

    [[nodiscard]] bool remember_frontier_state(
        const Task& task, WorkerContext& ctx, StateBestMap& best_gain) const
    {
        StateKey key = make_state_key(task, ctx);
        const auto it = best_gain.find(key);
        if (it != best_gain.end()) {
            if (it->second >= task.current_gain) return true;
            it->second = task.current_gain;
            return false;
        }
        best_gain.emplace(std::move(key), task.current_gain);
        return false;
    }

    [[nodiscard]] bool frontier_state_is_stale(
        const Task& task, WorkerContext& ctx, const StateBestMap& best_gain) const
    {
        StateKey key = make_state_key(task, ctx);
        const auto it = best_gain.find(key);
        return it != best_gain.end() && it->second > task.current_gain;
    }

    void dfs(Task& task, WorkerContext& ctx) {
        if (root_interval_closed()) return;
        node_count_.fetch_add(1, std::memory_order_relaxed);
        if (seen_search_state(task, ctx)) return;
        const long double ub = lagrangian_upper_bound(task, ctx);
        capture_root_cost_lb(ub);
        if (root_interval_closed()) return;
        if (gain_upper_bound_cannot_beat_best_cost(ub)) return;

        maybe_submit_greedy(task, ctx);
        if (root_interval_closed()) return;

        const int branch_line = choose_branch_line(task);
        if (branch_line < 0) {
            submit_candidate(task.current_gain, task.current_choice);
            return;
        }
        const std::size_t dual_mark = save_current_dual(ctx);
        const long double exclude_ub =
            exclude_branch_upper_bound(task, ctx, branch_line, ub);

        {
            const UndoMark undo = mark_undo(task, ctx);
            const int delta = positive_gain(task, branch_line);
            if (delta > 0) {
                select_line(task, ctx, branch_line);
                task.current_gain += delta;
                task.current_choice.push_back(branch_line);
                dfs(task, ctx);
                restore_saved_dual(ctx, dual_mark);
                task.current_choice.pop_back();
                task.current_gain -= delta;
            }
            undo_to(task, ctx, undo);
        }

        if (gain_upper_bound_cannot_beat_best_cost(exclude_ub))
            return;

        {
            const UndoMark undo = mark_undo(task, ctx);
            disable_line(task, ctx, branch_line);
            dfs(task, ctx);
            restore_saved_dual(ctx, dual_mark);
            undo_to(task, ctx, undo);
        }
    }

    void run_task(Task task, WorkerContext& ctx) {
        ensure_task_search_buffers(task);
        prepare_worker_context(ctx);
        if (root_interval_closed()) return;
        maybe_submit_greedy(task, ctx);
        if (root_interval_closed()) return;
        dfs(task, ctx);
    }

    void expand_task(const Task& task, std::vector<Task>& children) const {
        children.clear();
        const int branch_line = choose_branch_line(task);
        if (branch_line < 0) return;

        children.reserve(2);

        Task include = task;
        ensure_task_search_buffers(include);
        const int delta = positive_gain(include, branch_line);
        if (delta > 0) {
            include.current_gain += delta;
            include.current_choice.push_back(branch_line);
            select_line_without_undo(include, branch_line);
            children.push_back(std::move(include));
        }

        Task exclude = task;
        ensure_task_search_buffers(exclude);
        disable_line_without_undo(exclude, branch_line);
        children.push_back(std::move(exclude));
    }

    [[nodiscard]] std::vector<Task> build_frontier(Task root, WorkerContext& ctx) {
        std::vector<Task> frontier;
        frontier.reserve(std::max<std::size_t>(
            4, static_cast<std::size_t>(worker_count_) * 4U));
        frontier.push_back(std::move(root));

        std::vector<Task> children;
        const unsigned frontier_multiplier =
            current_best_cost() >= 93 ? 4U : 2U;
        const std::size_t target =
            std::max<std::size_t>(1, static_cast<std::size_t>(worker_count_) *
                                      static_cast<std::size_t>(frontier_multiplier));
        StateBestMap frontier_best;
        frontier_best.reserve(target * 4);
        static_cast<void>(remember_frontier_state(frontier.front(), ctx, frontier_best));

        while (frontier.size() < target) {
            if (root_interval_closed()) {
                frontier.clear();
                break;
            }
            const auto it = std::max_element(
                frontier.begin(), frontier.end(),
                [](const Task& a, const Task& b) {
                    if (a.active_count != b.active_count) return a.active_count < b.active_count;
                    return a.current_gain > b.current_gain;
                });
            if (it == frontier.end()) break;

            Task task = std::move(*it);
            frontier.erase(it);
            if (frontier_state_is_stale(task, ctx, frontier_best))
                continue;

            ctx.lagrangian.dual_ready = false;
            const long double ub = lagrangian_upper_bound(task, ctx);
            capture_root_cost_lb(ub);
            if (root_interval_closed()) {
                frontier.clear();
                break;
            }
            if (gain_upper_bound_cannot_beat_best_cost(ub)) {
                if (frontier.empty()) break;
                continue;
            }

            maybe_submit_greedy(task, ctx);
            if (root_interval_closed()) {
                frontier.clear();
                break;
            }

            const int branch_line = choose_branch_line(task);
            if (branch_line < 0) {
                submit_candidate(task.current_gain, task.current_choice);
                if (frontier.empty()) break;
                continue;
            }
            const long double exclude_ub =
                exclude_branch_upper_bound(task, ctx, branch_line, ub);

            children.clear();
            children.reserve(2);

            Task include = task;
            ensure_task_search_buffers(include);
            const int delta = positive_gain(include, branch_line);
            if (delta > 0) {
                include.current_gain += delta;
                include.current_choice.push_back(branch_line);
                select_line_without_undo(include, branch_line);
                children.push_back(std::move(include));
            }

            if (!gain_upper_bound_cannot_beat_best_cost(exclude_ub)) {
                Task exclude = task;
                ensure_task_search_buffers(exclude);
                disable_line_without_undo(exclude, branch_line);
                children.push_back(std::move(exclude));
            }

            for (Task& child : children) {
                if (remember_frontier_state(child, ctx, frontier_best))
                    continue;
                frontier.push_back(std::move(child));
            }
        }

        return frontier;
    }

    void solve_frontier_parallel(const std::vector<Task>& frontier) {
        const unsigned threads = std::min<unsigned>(
            worker_count_, static_cast<unsigned>(frontier.size()));
        if (threads <= 1) {
            WorkerContext ctx;
            run_task(frontier.front(), ctx);
            return;
        }

        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(threads - 1));

        auto worker = [&]() {
            WorkerContext ctx;
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
    }

    const std::vector<HeavyLine>&           lines_;
    const std::vector<std::vector<int>>&    incidence_;
    unsigned                                worker_count_ = 1;
    int                                     cur_n_ = 0;
    int                                     words_ = 0;
    std::vector<std::vector<int>>           ordered_point_lines_;
    std::size_t                             total_incidence_ = 0;
    std::atomic<int>                        best_gain_{0};
    std::atomic<bool>                       abort_{false};   // set by watchdog to cancel mid-run
    int                                     best_choice_gain_ = 0;
    int                                     root_cost_lb_ = 0;
    int                                     raw_root_cost_lb_ = 0;
    int                                     root_cost_lb_floor_ = 0;
    bool                                    root_lb_captured_ = false;
    std::atomic<long long>                  node_count_{0};
    std::vector<int>                        best_choice_;
    std::mutex                              best_choice_mutex_;
    static constexpr std::size_t            kSeenStateShardCount = 64;
    std::array<SeenStateShard, kSeenStateShardCount> seen_state_shards_{};
};

// ---------------------------------------------------------------------------
// Cover inequality lower bound
// Finds the maximum set S ⊆ [n] such that no heavy line covers more than 2
// points of S, then returns ⌈|S|/2⌉ as a valid lower bound on optimal cost.
//
// Proof of validity: every line in any cover (heavy, 2-point, or 1-point)
// covers at most 2 points of S by construction, so ⌈|S|/2⌉ lines are needed
// just to cover S.  Equivalently, this adds the set-cover inequality
//   ∑_{L : |L∩S| ≥ 1} x_L ≥ ⌈|S|/2⌉
// to the LP relaxation at the root, a tighter class of cut than the
// Lagrangian can certify.
//
// Algorithm: start with S = [n] and greedily remove the point that lies in
// the most currently-violated lines (lines with |L∩S| ≥ 3), updating
// violation counts incrementally until no line has more than 2 points in S.
// ---------------------------------------------------------------------------
[[nodiscard]] int cover_inequality_lb(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
{
    if (active_ids.empty() || n == 0) return 0;

    const int m = static_cast<int>(active_ids.size());

    // local_inc[p] = indices into active_ids for lines containing point p.
    std::vector<std::vector<int>> local_inc(static_cast<std::size_t>(n));
    std::vector<int> line_size(static_cast<std::size_t>(m), 0);
    for (int i = 0; i < m; ++i) {
        for (int p : lines[static_cast<std::size_t>(active_ids[i])].points) {
            if (p >= n) break;
            local_inc[static_cast<std::size_t>(p)].push_back(i);
            ++line_size[i];
        }
    }

    // S starts as all n points.  violation_count[p] = number of active lines
    // through p that currently have |L ∩ S| >= 3.
    std::vector<bool> in_S(static_cast<std::size_t>(n), true);
    int S_size = n;

    std::vector<int> violation_count(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < m; ++i) {
        if (line_size[i] >= 3) {
            for (int p : lines[static_cast<std::size_t>(active_ids[i])].points) {
                if (p >= n) break;
                ++violation_count[static_cast<std::size_t>(p)];
            }
        }
    }

    // Greedy: remove the point with the most violated-line incidences.
    while (true) {
        int worst_p = -1, worst_v = 0;
        for (int p = 0; p < n; ++p) {
            if (in_S[static_cast<std::size_t>(p)] &&
                violation_count[static_cast<std::size_t>(p)] > worst_v) {
                worst_v = violation_count[static_cast<std::size_t>(p)];
                worst_p = p;
            }
        }
        if (worst_v == 0) break;

        in_S[static_cast<std::size_t>(worst_p)] = false;
        --S_size;

        for (int i : local_inc[static_cast<std::size_t>(worst_p)]) {
            const bool was_violated = (line_size[i] >= 3);
            --line_size[i];
            if (was_violated && line_size[i] < 3) {
                // Line i just fell below the threshold; un-charge its S members.
                for (int q : lines[static_cast<std::size_t>(active_ids[i])].points) {
                    if (q >= n) break;
                    if (in_S[static_cast<std::size_t>(q)])
                        --violation_count[static_cast<std::size_t>(q)];
                }
            }
        }
    }

    return (S_size + 1) / 2;  // ceil(|S| / 2)
}

[[nodiscard]] int solve_exact(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines,
    const std::vector<std::vector<int>>& incidence, int root_lb_floor,
    const std::vector<int>& greedy_seed_line_ids, std::vector<int>& warm_start_solution,
    bool& timed_out, int& out_lb, int& out_lb_raw, long long& out_nodes)
{
    // Strengthen the root lower-bound floor with the cover inequality before
    // handing off to the gain solver.  If this closes the gap at the root,
    // root_interval_closed() fires after the first Lagrangian evaluation and
    // the DFS is eliminated entirely with zero branching.
    const int cover_lb = cover_inequality_lb(n, active_ids, lines);
    const int effective_lb_floor = std::max(root_lb_floor, cover_lb);
    ExactGainSolver solver(lines, incidence);
    const GainSolveResult result =
        solver.solve_for_n(
            n, active_ids, warm_start_solution, greedy_seed_line_ids, effective_lb_floor);
    warm_start_solution = result.chosen_line_ids;
    timed_out = result.timed_out;
    out_lb = result.root_cost_lb;
    out_lb_raw = result.raw_root_cost_lb;
    out_nodes = result.nodes;
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
    const auto total_start = std::chrono::steady_clock::now();
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
    long long total_nodes_processed = 0;
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
        const int ub_raw = greedy.cost;
        int incumbent = std::min(n, prev + 1);
        incumbent     = std::min(incumbent, greedy.cost);
        const int lb0_floor = warm_start_solution.empty() ? 0 : prev;
        const std::size_t active_count = active_ids.size();
        const int ub0_val = incumbent;

        bool timed_out = false;
        int out_lb = 0;
        int out_lb_raw = 0;
        long long out_nodes = 0;
        const int answer = solve_exact(
            n, active_ids, lines, incidence, lb0_floor, greedy.line_ids,
            warm_start_solution, timed_out, out_lb, out_lb_raw, out_nodes);
        const int gap = ub0_val - out_lb;
        const int gap_raw = ub_raw - out_lb_raw;

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();

        display.clear();                             // wipe status line, after timing
        total_nodes_processed += out_nodes;         // aggregate outside the timed region

        if (timed_out) {
            std::cout << "[stopped] N=" << n
                      << " time limit of " << kPerNTimeLimitSeconds
                      << "s exceeded (cancelled mid-run after "
                      << format_seconds(elapsed) << ")\n";
            break;
        }

        std::cout << "N="     << n
                  << " prime=" << primes[static_cast<std::size_t>(n - 1)]
                  << " active=" << active_count
                  << " ub0=" << ub0_val
                  << " lb=" << out_lb
                  << " gap=" << gap
                  << " ub_raw=" << ub_raw
                  << " lb_raw=" << out_lb_raw
                  << " gap_raw=" << gap_raw
                  << " lines=" << answer
                  << " time="  << format_seconds(elapsed)
                  << " nodes=" << out_nodes
                  << '\n';
        prev = answer;

        // Fallback post-hoc check: catches the (now unlikely) case where the
        // solve finished just before the watchdog fired but still ran over limit.
        if (kPerNTimeLimitSeconds > 0.0 && elapsed > kPerNTimeLimitSeconds) {
            std::cout << "[stopped] per-N time limit of " << kPerNTimeLimitSeconds
                      << "s exceeded\n";
            break;
        }
    }
    const double total_elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - total_start).count();
    std::cout << "total time=" << format_seconds(total_elapsed) << " total nodes=" << total_nodes_processed << '\n';
    std::cout << '\n'
              << "Output fields:\n"
              << "  N=      — The integer being solved for.\n"
              << "  prime=  — The value of the N-th prime number.\n"
              << "  active= — The number of heavy-line candidates (lines through 3 or more points) entering the exact search.\n"
              << "  ub0=    — The best-known solution cost at the start of the search, before any branching.\n"
              << "  lb=     — The Lagrangian cost lower bound at the root node, before any branching.\n"
              << "  gap=    — ub0 minus lb; the proof burden the search must close to certify optimality.\n"
              << "  ub_raw= — The greedy upper bound before the monotonic cap (min with prev_opt+1) is applied.\n"
              << "  lb_raw= — The Lagrangian cost lower bound at the root before the monotonic floor (max with prev_opt) is applied.\n"
              << "  gap_raw=— ub_raw minus lb_raw; the raw proof burden before tightening. Shows intrinsic algorithm performance.\n"
              << "  lines=  — The proven optimal number of lines covering all N points — the answer.\n"
              << "  time=   — Wall-clock seconds elapsed for the greedy construction and full exact search for this N.\n"
              << "  nodes=  — Total calls to the DFS function; a direct measure of how hard the instance was to solve.\n";
    return 0;
}
