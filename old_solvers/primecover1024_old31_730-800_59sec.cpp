/*
 * primecover1024.cpp
 *
 * Exact minimum line cover solver for prime points (i, p_i).
 *
 * Live solver:
 * - Enumerates every heavy line containing at least 3 points.
 * - Solves the residual exact-gain search with 1024-bit masks, greedy
 * seeding, a Lagrangian upper bound, and a shared atomic incumbent.
 * - Applies a root cover-inequality lower bound before branch-and-bound.
 * - Closes any uncovered residue with exact 2-point / 1-point lines at
 * cost ceil(residual / 2).
 *  
 * Experiment Log:
 *   1. We tested an early-stop heuristic in the Lagrangian subgradient loop that was meant to cut bound time by exiting after stretches of non-improving iterations, with variants at 8 steps, 16 steps, and a delayed 24/12 trigger.
 *      In measured runs from N=730..800, every variant regressed against the kept specialized-eval + adaptive-iteration build: totals rose from about 60.19s/60.20s to 81.99s, 63.47s, and 60.96s, and node counts also jumped sharply.
 *      The simplest reason is that the earlier exits weakened the bound enough to expand far more search nodes than the saved dual iterations could repay.
 *      Future attempts should avoid blanket early termination and instead use stronger convergence signals or tighter gating that preserves heavy-tail bound quality.
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <print>
#include <ranges>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std::chrono_literals;

namespace {

std::atomic<int>* forced_atomic_ptr = nullptr;
std::atomic<long long>* lag_iters_atomic_ptr = nullptr;
std::atomic<int>* depth_max_atomic_ptr = nullptr;
std::atomic<int>* seen_prune_atomic_ptr = nullptr;
std::atomic<int>* lag_prune_atomic_ptr = nullptr;
std::atomic<int>* strong_branch_atomic_ptr = nullptr;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
namespace config {
    constexpr int    kBitCapacity          = 1024;
    constexpr int    kBitWords             = kBitCapacity / 64;
    constexpr int    kStartN               = 730;   // N to start from in the timed sweep; set to 0 to run from 1
    constexpr int    kExecutionLimit       = 800;   // hard cap and default run-to N
    // Maximum wall-clock seconds allowed per N. Set to 0.0 to disable (no limit). 60=1 min
    constexpr double kPerNTimeLimitSeconds = 25;
} // namespace config

// ---------------------------------------------------------------------------
// Utility Helpers
// ---------------------------------------------------------------------------
namespace util {
    template <typename T>
    [[nodiscard]] constexpr std::size_t as_size(T v) noexcept {
        return static_cast<std::size_t>(v);
    }

    template <typename T>
    [[nodiscard]] constexpr int as_int(T v) noexcept {
        return static_cast<int>(v);
    }

    [[nodiscard]] int words_for_n(int n) noexcept { return (n + 63) / 64; }
    [[nodiscard]] int ceil_half(int v)   noexcept { return (v + 1) / 2; }

    [[nodiscard]] std::string format_seconds(double s) {
        return std::format("{:.{}f}s", s, s < 10.0 ? 6 : 3);
    }

    [[nodiscard]] std::uint64_t mix_u64(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    [[nodiscard]] std::int64_t abs64(std::int64_t x) noexcept { return x < 0 ? -x : x; }
} // namespace util

// ---------------------------------------------------------------------------
// Debug & Tracing
// ---------------------------------------------------------------------------
namespace debug {
    struct DebugOptions {
        int         trace_n = 0;
        bool        force_exact = false;
        std::string trace_prefix;
    };

    [[nodiscard]] const DebugOptions& debug_options() {
        static const DebugOptions options = [] {
            DebugOptions out;
            if (const char* value = std::getenv("PRIMECOVER_TRACE_N"))
                out.trace_n = std::atoi(value);
            if (const char* value = std::getenv("PRIMECOVER_FORCE_EXACT"))
                out.force_exact = std::atoi(value) != 0;
            if (const char* value = std::getenv("PRIMECOVER_TRACE_PREFIX"))
                out.trace_prefix = value;
            return out;
        }();
        return options;
    }

    [[nodiscard]] bool debug_trace_enabled_for_n(int n) {
        const auto& options = debug_options();
        return options.trace_n > 0 && options.trace_n == n && !options.trace_prefix.empty();
    }

    [[nodiscard]] std::string debug_trace_path(const char* suffix) {
        return debug_options().trace_prefix + suffix;
    }

    void debug_write_text(const char* suffix, const std::string& text, const char* mode) {
        const auto& options = debug_options();
        if (options.trace_prefix.empty()) return;
        if (std::FILE* file = std::fopen(debug_trace_path(suffix).c_str(), mode)) {
            std::fwrite(text.data(), 1, text.size(), file);
            std::fclose(file);
        }
    }

    void debug_truncate_text(const char* suffix) {
        debug_write_text(suffix, "", "w");
    }

    void debug_append_text(const char* suffix, const std::string& text) {
        static std::mutex mutex;
        std::scoped_lock lock(mutex);
        debug_write_text(suffix, text, "a");
    }
} // namespace debug

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------
namespace geom {
    using util::as_size;
    using util::mix_u64;
    using util::abs64;

    struct LineKey {
        std::int64_t a{}, b{}, c{};
        bool operator==(const LineKey&) const = default;
    };

    struct SlopeKey {
        int dy = 0;
        int dx = 1;
        bool operator==(const SlopeKey&) const = default;
    };

    struct SlopeKeyHash {
        std::size_t operator()(const SlopeKey& k) const noexcept {
            const std::uint64_t packed =
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(k.dy)) << 32) |
                static_cast<std::uint32_t>(k.dx);
            return as_size(mix_u64(packed));
        }
    };

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

    [[nodiscard]] SlopeKey reduced_slope(int dx, int dy) noexcept {
        const int g = std::gcd(dx, dy);
        return {dy / g, dx / g};
    }

    [[nodiscard]] bool line_contains_point(const LineKey& line, int x, int y) noexcept {
        return line.a * static_cast<std::int64_t>(x) +
               line.b * static_cast<std::int64_t>(y) + line.c == 0;
    }
} // namespace geom

// ---------------------------------------------------------------------------
// Core Data & Generation
// ---------------------------------------------------------------------------
namespace core {
    using namespace geom;
    using util::as_size;
    using util::as_int;
    using util::ceil_half;
    using util::words_for_n;

    struct BitMask1024 {
        std::array<std::uint64_t, config::kBitWords> words{};

        void set(int bit) noexcept   { words[bit >> 6] |=  (1ULL << (bit & 63)); }
        void reset(int bit) noexcept { words[bit >> 6] &= ~(1ULL << (bit & 63)); }
        bool test(int bit) const noexcept {
            return ((words[bit >> 6] >> (bit & 63)) & 1ULL) != 0;
        }

        int intersection_count(const BitMask1024& other, int used_words) const noexcept {
            int total = 0;
            for (int i = 0; i < used_words; ++i)
                total += std::popcount(words[i] & other.words[i]);
            return total;
        }

        BitMask1024 and_not(const BitMask1024& other) const noexcept {
            BitMask1024 out;
            for (int i = 0; i < config::kBitWords; ++i) out.words[i] = words[i] & ~other.words[i];
            return out;
        }

        static BitMask1024 prefix(int n) noexcept {
            BitMask1024 out;
            const int full_words = n / 64;
            for (int i = 0; i < full_words; ++i)
                out.words[i] = std::numeric_limits<std::uint64_t>::max();
            if (const int tail = n % 64; tail != 0 && full_words < config::kBitWords)
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
    };

    struct HeavyLine {
        BitMask1024      mask;
        std::vector<int> points;
        int              activate_at = 0;
    };

    struct GreedyCoverSolution {
        int              cost = 0;
        std::vector<int> line_ids;
    };

    struct GainSolveResult {
        int              best_gain = 0;
        int              root_cost_lb = 0;
        int              raw_root_cost_lb = 0;
        int              root_productive_lines = 0;      // productive heavy lines at the post-reduction root
        int              root_productive_incidence = 0;  // sum of root line_cover[id] over those productive lines
        long long        nodes = 0;
        std::vector<int> chosen_line_ids;
        bool             timed_out = false;
    };

    struct ExactSolveResult {
        int              answer = 0;
        int              lb = 0;
        int              lb_raw = 0;
        int              lb_cov = 0;     // root lower bound from cover_inequality_lb()
        int              lb_lp = 0;      // reserved output field; no LP lower bound is active
        int              prod = 0;       // productive heavy lines in the post-reduction root task
        int              pinc = 0;       // sum of line_cover over those productive root lines
        long long        nodes = 0;
        std::vector<int> chosen_line_ids;
        bool             timed_out = false;
    };

    struct WitnessCover {
        std::vector<LineKey> lines;
    };

    struct ChoiceReplayStep {
        int depth = 0;
        int line_id = -1;
        int first_point = -1;
        int second_point = -1;
        int point_count = 0;
        int marginal = 0;
        int cumulative = 0;
    };

    struct ChoiceReplay {
        std::vector<ChoiceReplayStep> steps;
        std::vector<int>              residual_points;
        int                           covered_points = 0;
        int                           residual_lines = 0;
    };

    [[nodiscard]] std::vector<int> generate_primes(int count) {
        if (count <= 0) return {};
        const double n = static_cast<double>(count);
        const int limit = (count < 10)
            ? 64
            : as_int(n * (std::log(n) + std::log(std::log(n)) + 3.0)) + 256;
        std::vector<bool> is_p(as_size(limit + 1), true);
        is_p[0] = is_p[1] = false;
        for (int p = 2; p * p <= limit; ++p)
            if (is_p[as_size(p)])
                for (int q = p * p; q <= limit; q += p)
                    is_p[as_size(q)] = false;
        std::vector<int> primes;
        primes.reserve(as_size(count));
        for (int v = 2; v <= limit && as_int(primes.size()) < count; ++v)
            if (is_p[as_size(v)]) primes.push_back(v);
        return primes;
    }

    [[nodiscard]] bool anchor_is_leftmost_for_slope(
        const std::vector<int>& primes, int anchor, const SlopeKey& slope) noexcept
    {
        const int step_x = slope.dx;
        const std::int64_t step_y = slope.dy;
        std::int64_t expected_y = static_cast<std::int64_t>(primes[as_size(anchor)]) - step_y;
        for (int x = anchor - step_x; x >= 0; x -= step_x, expected_y -= step_y) {
            if (expected_y == static_cast<std::int64_t>(primes[as_size(x)])) return false;
        }
        return true;
    }

    [[nodiscard]] std::vector<HeavyLine> enumerate_heavy_lines(const std::vector<int>& primes) {
        const int n = as_int(primes.size());
        std::vector<HeavyLine> lines;
        lines.reserve(as_size(n));

        std::unordered_map<SlopeKey, std::vector<int>, SlopeKeyHash> slope_buckets;
        slope_buckets.reserve(as_size(n));
        for (int i = 0; i < n; ++i) {
            slope_buckets.clear();
            const int base_prime = primes[as_size(i)];
            for (int j = i + 1; j < n; ++j) {
                const SlopeKey slope = reduced_slope(j - i, primes[as_size(j)] - base_prime);
                slope_buckets[slope].push_back(j);
            }

            for (auto& [slope, bucket] : slope_buckets) {
                if (bucket.size() < 2) continue;
                if (!anchor_is_leftmost_for_slope(primes, i, slope)) continue;

                HeavyLine hl;
                hl.points.reserve(bucket.size() + 1);
                hl.points.push_back(i);
                hl.mask.set(i);
                for (int p : bucket) {
                    hl.mask.set(p);
                    hl.points.push_back(p);
                }
                hl.activate_at = hl.points[2] + 1;
                lines.push_back(std::move(hl));
            }
        }
        std::ranges::sort(lines, [](const HeavyLine& a, const HeavyLine& b) {
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
                const int c = lines[as_size(id)].mask.intersection_count(active, words);
                if (c > bc || (c == bc && best != -1 && lines[as_size(id)].points.size() > lines[as_size(best)].points.size())) { 
                    bc = c; best = id; 
                }
            }
            if (best < 0) break;
            active = active.and_not(lines[as_size(best)].mask);
            count -= bc;
            ++result.cost;
            result.line_ids.push_back(best);
        }
        result.cost += ceil_half(count);
        return result;
    }

    [[nodiscard]] int cost_from_gain(int n, int gain) noexcept {
        return ceil_half(n - gain);
    }

    [[nodiscard]] WitnessCover build_witness_cover(
        int n, const std::vector<int>& primes, const std::vector<HeavyLine>& lines, const std::vector<int>& chosen_heavy_line_ids)
    {
        WitnessCover witness;
        witness.lines.reserve(chosen_heavy_line_ids.size() + as_size(ceil_half(n)));

        BitMask1024 covered;
        for (int id : chosen_heavy_line_ids) {
            const auto& heavy = lines[as_size(id)];
            witness.lines.push_back(canonical_line(
                heavy.points[0], primes[as_size(heavy.points[0])],
                heavy.points[1], primes[as_size(heavy.points[1])]));
            for (int p : heavy.points) {
                if (p >= n) break;
                covered.set(p);
            }
        }

        int pending = -1;
        for (int p = 0; p < n; ++p) {
            if (covered.test(p)) continue;
            if (pending < 0) {
                pending = p;
                continue;
            }
            witness.lines.push_back(canonical_line(
                pending, primes[as_size(pending)],
                p,       primes[as_size(p)]));
            pending = -1;
        }
        if (pending >= 0)
            witness.lines.push_back({1, 0, -static_cast<std::int64_t>(pending)});

        return witness;
    }

    [[nodiscard]] bool witness_covers_point(const WitnessCover& witness, int x, int y) noexcept {
        for (const LineKey& line : witness.lines) {
            if (line_contains_point(line, x, y)) return true;
        }
        return false;
    }

    [[nodiscard]] ChoiceReplay replay_choice_trace(
        int n, const std::vector<HeavyLine>& lines, const std::vector<int>& chosen_line_ids)
    {
        ChoiceReplay replay;
        replay.steps.reserve(chosen_line_ids.size());

        BitMask1024 active = BitMask1024::prefix(n);
        int covered = 0;
        for (std::size_t depth = 0; depth < chosen_line_ids.size(); ++depth) {
            const int id = chosen_line_ids[depth];
            const auto& line = lines[as_size(id)];

            int marginal = 0;
            for (int p : line.points) {
                if (p >= n) break;
                if (!active.test(p)) continue;
                active.reset(p);
                ++marginal;
            }
            covered += marginal;
            replay.steps.push_back({
                as_int(depth + 1), id,
                line.points.empty() ? -1 : line.points[0],
                line.points.size() > 1 ? line.points[1] : -1,
                as_int(line.points.size()), marginal, covered
            });
        }

        active.for_each_set_bit(words_for_n(n), [&](int p) { replay.residual_points.push_back(p); });
        replay.covered_points = covered;
        replay.residual_lines = ceil_half(as_int(replay.residual_points.size()));
        return replay;
    }

    [[nodiscard]] std::string format_choice_trace(
        int n, int answer, const ChoiceReplay& replay, const std::vector<int>& chosen_line_ids)
    {
        std::string out = std::format("n={} answer={} heavy_lines={} covered={} residual={} residual_lines={}\n",
            n, answer, chosen_line_ids.size(), replay.covered_points, replay.residual_points.size(), replay.residual_lines);
        for (const auto& step : replay.steps) {
            out += std::format("depth={} line={} first={} second={} size={} marginal={} cumulative={}\n",
                step.depth, step.line_id, step.first_point, step.second_point, step.point_count, step.marginal, step.cumulative);
        }
        out += "residual_points=";
        for (std::size_t i = 0; i < replay.residual_points.size(); ++i) {
            out += std::format("{}{}", i != 0 ? "," : "", replay.residual_points[i]);
        }
        out += '\n';
        return out;
    }
} // namespace core

// ---------------------------------------------------------------------------
// Solver Engine
// ---------------------------------------------------------------------------
namespace solver {
    using namespace core;
    using util::as_size;
    using util::as_int;
    using util::words_for_n;
    using util::mix_u64;
    using debug::debug_options;
    using debug::debug_trace_enabled_for_n;
    using debug::debug_truncate_text;
    using debug::debug_append_text;

class ExactGainSolver {
public:
    ExactGainSolver(const std::vector<HeavyLine>& lines, const std::vector<std::vector<int>>& incidence)
        : lines_(lines), incidence_(incidence)
        , worker_count_(std::max(1U, std::thread::hardware_concurrency()))
        , total_incidence_(std::ranges::fold_left(lines, std::size_t{0},
              [](std::size_t total, const HeavyLine& line) { return total + line.points.size(); })) {}

    [[nodiscard]] GainSolveResult solve_for_n(
        int n, const std::vector<int>& active_ids, const std::vector<int>& warm_start_line_ids,
        const std::vector<int>& greedy_seed_line_ids, int root_cost_lb_floor)
    {
        cur_n_  = n;
        words_  = words_for_n(n);
        best_gain_.store(0, std::memory_order_relaxed);
        abort_.store(false, std::memory_order_relaxed);
        trace_this_n_ = debug_trace_enabled_for_n(n);
        force_exact_forced_lines_ = trace_this_n_ && debug_options().force_exact;
        
        if (trace_this_n_) {
            debug_truncate_text(".events.txt");
            debug_truncate_text(".chosen.txt");
        }
        
        root_cost_lb_ = 0;
        raw_root_cost_lb_ = 0;
        root_cost_lb_floor_ = root_cost_lb_floor;
        root_lb_captured_ = false;
        node_count_.store(0, std::memory_order_relaxed);
        best_choice_gain_ = 0;
        best_choice_.clear();
        best_choice_.reserve(lines_.size());
        clear_seen_state_table();

        ordered_point_lines_.assign(as_size(n), {});
        for (int p = 0; p < n; ++p) {
            auto& ordered = ordered_point_lines_[as_size(p)];
            ordered.reserve(incidence_[as_size(p)].size());
            for (int id : incidence_[as_size(p)]) {
                if (lines_[as_size(id)].activate_at <= n) ordered.push_back(id);
            }
            std::ranges::stable_sort(ordered, [&](int lhs, int rhs) {
                const auto sa = lines_[as_size(lhs)].points.size();
                const auto sb = lines_[as_size(rhs)].points.size();
                if (sa != sb) return sa > sb;
                return lhs < rhs;
            });
        }

        Task root = make_root_task(n, active_ids);
        WorkerContext root_ctx;
        prepare_worker_context(root_ctx);
        seed_from_line_ids(root, warm_start_line_ids);
        seed_from_line_ids(root, greedy_seed_line_ids);
        apply_forced_line_reduction(root, root_ctx, false);
        
        const auto [root_productive_lines, root_productive_incidence] = root_productive_metrics(root);
        const int greedy_gain = greedy_completion(root, root_ctx, true);
        if (greedy_gain > current_best_gain()) submit_candidate(greedy_gain, root_ctx.greedy_lines);
        greedy_complete_from_seed(root, root_ctx, warm_start_line_ids);

        // Watchdog thread: sets abort_ once the per-N wall-clock limit is reached.
        // The flag is consumed inside lagrangian_upper_bound() — one relaxed load
        // amortised over all the expensive bound work — so there is zero overhead on
        // the DFS branch/select/undo hot path during normal operation.
        std::atomic<bool> watchdog_done{false};
        std::thread watchdog_thread;
        if (config::kPerNTimeLimitSeconds > 0.0) {
            watchdog_thread = std::thread([&]() {
                const auto deadline = std::chrono::steady_clock::now() + std::chrono::duration<double>(config::kPerNTimeLimitSeconds);
                while (!watchdog_done.load(std::memory_order_relaxed)) {
                    const auto now = std::chrono::steady_clock::now();
                    if (now >= deadline) {
                        abort_.store(true, std::memory_order_relaxed);
                        return;
                    }
                    const auto remaining = deadline - now;
                    const auto slice = std::chrono::duration_cast<std::chrono::steady_clock::duration>(50ms);
                    std::this_thread::sleep_for(remaining < slice ? remaining : slice);
                }
            });
        }

        if (worker_count_ <= 1) {
            run_task(std::move(root), root_ctx);
        } else {
            std::vector<Task> frontier = build_frontier(std::move(root), root_ctx);
            if (frontier.size() == 1) {
                run_task(std::move(frontier.front()), root_ctx);
            } else if (!frontier.empty()) {
                solve_frontier_parallel(frontier);
            }
        }

        watchdog_done.store(true, std::memory_order_relaxed);
        if (watchdog_thread.joinable()) watchdog_thread.join();

        GainSolveResult out;
        out.best_gain = current_best_gain();
        out.root_cost_lb = root_cost_lb_;
        out.raw_root_cost_lb = raw_root_cost_lb_;
        out.root_productive_lines = root_productive_lines;
        out.root_productive_incidence = root_productive_incidence;
        out.nodes = node_count_.load(std::memory_order_relaxed);
        out.chosen_line_ids = best_choice_;
        out.timed_out = abort_.load(std::memory_order_relaxed);
        return out;
    }

private:
    using HeapItem = std::pair<int, int>;
    using DualValue = double;
    using CompactPointIndex = std::uint16_t;

    static_assert(config::kBitCapacity <= std::numeric_limits<CompactPointIndex>::max());

    struct CompactLine3 {
        CompactPointIndex a, b, c;
    };

    struct CompactLine4 {
        CompactPointIndex a, b, c, d;
    };

    struct CompactLine5 {
        CompactPointIndex a, b, c, d, e;
    };

    struct Task {
        BitMask1024                   active;
        int                           active_count  = 0;
        int                           current_gain  = 0;
        std::vector<std::uint16_t>    line_cover;
        std::vector<int>              productive_degree;
        std::vector<int>              sole_productive_line;
        std::vector<unsigned char>    dirty_points;
        std::vector<unsigned char>    available;
        std::vector<int>              current_choice;
        // Compact set of productive line IDs: available[id]==1 && line_cover[id]>=3.
        // productive_pos[id] = index in productive_ids, or -1 if not in set.
        std::vector<int>              productive_ids;
        std::vector<int>              productive_pos;
        // Compact set of unavailable lines whose current residual coverage is
        // still gain-positive.  Duplicate-state pruning keys only on this set:
        // once a blocked line falls below coverage 3 it can never become
        // relevant again because the active set only shrinks.
        std::vector<int>              blocked_productive_ids;
        std::vector<int>              blocked_productive_pos;
    };

    struct UndoMark {
        std::size_t disabled_lines = 0;
        std::size_t decremented_lines = 0;
        std::size_t newly_covered_points = 0;
        std::size_t covered_point_metadata = 0;
        std::size_t removed_from_productive = 0;
        std::size_t removed_from_blocked_productive = 0;
    };

    struct StateKey {
        BitMask1024      active;
        std::vector<int> blocked_productive_ids;

        bool operator==(const StateKey& other) const noexcept {
            return active.words == other.active.words && blocked_productive_ids == other.blocked_productive_ids;
        }
    };

    struct StateKeyHash {
        std::size_t operator()(const StateKey& key) const noexcept {
            std::uint64_t h = 0x9e3779b97f4a7c15ULL;
            for (std::uint64_t word : key.active.words)
                h ^= mix_u64(word + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            for (int id : key.blocked_productive_ids)
                h ^= mix_u64(static_cast<std::uint64_t>(id) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            return as_size(h);
        }
    };

    using StateBestMap = std::unordered_map<StateKey, int, StateKeyHash>;

    struct SeenStateShard {
        std::mutex   mutex;
        StateBestMap best_gain;
    };

    struct LagrangianScratch {
        std::vector<int>         active_points;
        std::vector<DualValue>   active_y;
        std::vector<DualValue>   active_best_y;
        std::vector<DualValue>   best_y;
        std::vector<DualValue>   active_subgrad;
        std::vector<DualValue>   work_line_slack;
        std::vector<long double> line_slack;
        std::vector<int>         compact_line_ids;
        std::vector<std::size_t> compact_line_offsets;
        std::vector<CompactPointIndex> compact_line_points;
        std::vector<CompactPointIndex> compact_point_aps;
        std::vector<CompactLine3> compact_eval3;
        std::vector<CompactLine4> compact_eval4;
        std::vector<CompactLine5> compact_eval5;
        std::vector<int>         compact_eval_other_lines;
        bool                     dual_ready = false;

        // Compact inverse: active-point index (ap) → productive-line compact indices (li).
        // Built alongside the forward compact in build_compact_productive_lines.
        // Eliminates the global-incidence + filter scan in the coordinate polish.
        std::vector<int>         compact_pt_prod_offsets;  // size: active_points.size()+1
        std::vector<int>         compact_pt_prod_flat;     // li values, total size = pinc
        std::vector<int>         compact_ap_fill_pos;      // scratch write cursors during build
        std::vector<int>         point_to_ap_scratch;      // global point → ap index; stays -1 except during build

        // Compact-indexed slack array: avoids random global-id slack lookups in the
        // polish hot path and lets the compiler vectorise more aggressively.
        std::vector<DualValue>   work_line_slack_compact;  // indexed by li (compact line index)

        // Precomputed (line_cover[id] - 2) for each compact productive line.
        // line_cover is a ~16 KB uint16_t array accessed with random ids; keeping a
        // compact sequential copy here eliminates that random-access load from the
        // eval_dense_fast/eval_dense_exact/polish-init hot loops and shrinks the
        // per-bound L1 working set from ~34 KB to ~18 KB.
        std::vector<int>         compact_rem2;             // indexed by li (compact line index)

        void ensure(std::size_t line_count, int current_n) {
            const std::size_t n = as_size(current_n);
            if (active_points.capacity() < n) active_points.reserve(n);
            if (active_y.size() != n) active_y.resize(n);
            if (active_best_y.size() != n) active_best_y.resize(n);
            if (best_y.size() != n) best_y.resize(n);
            if (active_subgrad.size() != n) active_subgrad.resize(n);
            if (work_line_slack.size() != line_count) work_line_slack.resize(line_count);
            if (line_slack.size() != line_count) line_slack.resize(line_count);
            if (compact_line_ids.capacity() < line_count) compact_line_ids.reserve(line_count);
            if (compact_line_offsets.capacity() < line_count + 1) compact_line_offsets.reserve(line_count + 1);
            if (compact_eval3.capacity() < line_count) compact_eval3.reserve(line_count);
            if (compact_eval4.capacity() < line_count) compact_eval4.reserve(line_count);
            if (compact_eval5.capacity() < line_count) compact_eval5.reserve(line_count);
            if (compact_eval_other_lines.capacity() < line_count) compact_eval_other_lines.reserve(line_count);
            // Ensure point_to_ap_scratch is always all-(-1) so the build helper
            // can write into it and cleanly reset just the touched entries.
            if (as_int(point_to_ap_scratch.size()) < current_n) point_to_ap_scratch.assign(as_size(current_n), -1);
        }
    };

    struct WorkerContext {
        LagrangianScratch          lagrangian;
        std::vector<HeapItem>      greedy_heap;
        std::array<int, config::kBitCapacity + 1> greedy_bucket_head{};
        std::vector<int>           greedy_bucket_next;
        std::vector<int>           greedy_lines;
        std::vector<int>           candidate_choice;
        std::vector<int>           disabled_lines;
        std::vector<int>           decremented_lines;
        std::vector<int>           newly_covered_points;
        std::vector<int>           removed_from_productive;
        std::vector<int>           removed_from_blocked_productive;
        std::vector<int>           state_blocked_ids;
        std::vector<int>           saved_dual_points;
        std::vector<DualValue>     saved_dual_values;
        std::vector<int>           forced_unique_count;
        std::vector<int>           forced_dirty_unique_count;
        std::vector<int>           forced_touched_lines;
        std::vector<int>           covered_point_degree;
        std::vector<int>           covered_point_sole;
        std::vector<unsigned char> covered_point_dirty;
    };

    static constexpr int kNoProductiveLine = -1;
    static constexpr int kManyProductiveLines = -2;
    static constexpr int kStrongBranchShortlist = 4;
    static constexpr int kStrongBranchMaxIters = 24;
    static constexpr int kStrongBranchMaxDepth = 6;

    struct BranchChoice {
        int         line = -1;
        long double exclude_ub = 0.0L;
    };

    struct ForcedChoice {
        int best_line = -1;
        int best_unique = 0;
        int best_gain = -1;
    };

    [[nodiscard]] Task make_root_task(int n, const std::vector<int>& active_ids) const {
        Task root;
        root.active = BitMask1024::prefix(n);
        root.active_count = n;
        root.line_cover.assign(lines_.size(), 0);
        root.productive_degree.assign(as_size(n), 0);
        root.sole_productive_line.assign(as_size(n), kNoProductiveLine);
        root.dirty_points.assign(as_size(n), 0);
        root.available.assign(lines_.size(), 0);
        root.productive_pos.assign(lines_.size(), -1);
        root.blocked_productive_pos.assign(lines_.size(), -1);
        root.current_choice.reserve(lines_.size());
        root.productive_ids.reserve(active_ids.size());
        root.blocked_productive_ids.reserve(64);
        
        for (int id : active_ids) {
            const int cover = lines_[as_size(id)].mask.intersection_count(root.active, words_);
            root.line_cover[as_size(id)] = static_cast<std::uint16_t>(cover);
            root.available[as_size(id)] = 1;
            if (cover >= 3) {
                root.productive_pos[as_size(id)] = as_int(root.productive_ids.size());
                root.productive_ids.push_back(id);
            }
        }
        rebuild_productive_point_metadata(root);
        return root;
    }

    void ensure_task_search_buffers(Task& task) const {
        if (task.current_choice.capacity() < lines_.size()) task.current_choice.reserve(lines_.size());
        if (task.blocked_productive_ids.capacity() < 64) task.blocked_productive_ids.reserve(64);
    }

    void prepare_worker_context(WorkerContext& ctx) const {
        ctx.lagrangian.ensure(lines_.size(), cur_n_);
        ctx.greedy_heap.clear();        ctx.greedy_heap.reserve(lines_.size());
        if (ctx.greedy_bucket_next.size() != lines_.size()) ctx.greedy_bucket_next.resize(lines_.size());
        ctx.greedy_lines.clear();       ctx.greedy_lines.reserve(lines_.size());
        ctx.candidate_choice.clear();   ctx.candidate_choice.reserve(lines_.size());
        ctx.disabled_lines.clear();     ctx.disabled_lines.reserve(lines_.size());
        ctx.decremented_lines.clear();  ctx.decremented_lines.reserve(total_incidence_);
        ctx.newly_covered_points.clear(); ctx.newly_covered_points.reserve(as_size(cur_n_));
        ctx.removed_from_productive.clear(); ctx.removed_from_productive.reserve(lines_.size());
        ctx.removed_from_blocked_productive.clear();
        ctx.removed_from_blocked_productive.reserve(lines_.size());
        ctx.state_blocked_ids.clear();  ctx.state_blocked_ids.reserve(64);
        ctx.saved_dual_points.clear();  ctx.saved_dual_points.reserve(as_size(cur_n_));
        ctx.saved_dual_values.clear();  ctx.saved_dual_values.reserve(as_size(cur_n_));
        
        if (ctx.forced_unique_count.size() != lines_.size()) ctx.forced_unique_count.resize(lines_.size());
        if (ctx.forced_dirty_unique_count.size() != lines_.size()) ctx.forced_dirty_unique_count.resize(lines_.size());
        ctx.forced_touched_lines.clear(); ctx.forced_touched_lines.reserve(lines_.size());
        ctx.covered_point_degree.clear(); ctx.covered_point_degree.reserve(as_size(cur_n_));
        ctx.covered_point_sole.clear(); ctx.covered_point_sole.reserve(as_size(cur_n_));
        ctx.covered_point_dirty.clear(); ctx.covered_point_dirty.reserve(as_size(cur_n_));
        ctx.lagrangian.dual_ready = false;
    }

    [[nodiscard]] int current_best_gain() const noexcept { return best_gain_.load(std::memory_order_relaxed); }
    [[nodiscard]] int current_best_cost() const noexcept { return cost_from_gain(cur_n_, current_best_gain()); }
    [[nodiscard]] bool root_interval_closed() const noexcept { return root_lb_captured_ && current_best_cost() <= root_cost_lb_; }
    [[nodiscard]] int gain_needed_to_beat_best_cost() const noexcept { return cur_n_ - 2 * current_best_cost() + 2; }
    [[nodiscard]] bool gain_upper_bound_cannot_beat_best_cost(long double total_gain_ub) const noexcept {
        return as_int(std::floor(total_gain_ub + 1e-9L)) < gain_needed_to_beat_best_cost();
    }

    void submit_candidate(int gain, const std::vector<int>& choice) {
        int observed = current_best_gain();
        while (gain > observed && !best_gain_.compare_exchange_weak(observed, gain, std::memory_order_relaxed, std::memory_order_relaxed)) {}
        if (gain < current_best_gain()) return;
        std::scoped_lock lock(best_choice_mutex_);
        if (gain > best_choice_gain_) {
            best_choice_gain_ = gain;
            best_choice_ = choice;
        }
    }

    [[nodiscard]] int positive_gain(const Task& task, int id) const noexcept {
        if (!task.available[as_size(id)]) return 0;
        const int gain = as_int(task.line_cover[as_size(id)]) - 2;
        return gain > 0 ? gain : 0;
    }

    [[nodiscard]] static UndoMark mark_undo(const WorkerContext& ctx) noexcept {
        return { ctx.disabled_lines.size(), ctx.decremented_lines.size(), ctx.newly_covered_points.size(),
                 ctx.covered_point_degree.size(), ctx.removed_from_productive.size(),
                 ctx.removed_from_blocked_productive.size() };
    }

    [[nodiscard]] std::size_t save_dual_for_points(WorkerContext& ctx, const std::vector<int>& points) const {
        const std::size_t mark = ctx.saved_dual_points.size();
        if (!ctx.lagrangian.dual_ready) return mark;
        for (int p : points) {
            ctx.saved_dual_points.push_back(p);
            ctx.saved_dual_values.push_back(ctx.lagrangian.best_y[as_size(p)]);
        }
        return mark;
    }

    [[nodiscard]] std::size_t save_current_dual(WorkerContext& ctx) const { return save_dual_for_points(ctx, ctx.lagrangian.active_points); }

    void restore_saved_dual(WorkerContext& ctx, std::size_t mark) const {
        for (std::size_t i = mark; i < ctx.saved_dual_points.size(); ++i) {
            ctx.lagrangian.best_y[as_size(ctx.saved_dual_points[i])] = ctx.saved_dual_values[i];
        }
        ctx.saved_dual_points.resize(mark);
        ctx.saved_dual_values.resize(mark);
        ctx.lagrangian.dual_ready = true;
    }

    void rebuild_productive_point_metadata(Task& task) const {
        std::ranges::fill(task.productive_degree, 0);
        std::ranges::fill(task.sole_productive_line, kNoProductiveLine);
        for (int id : task.productive_ids) {
            for (int p : lines_[as_size(id)].points) {
                if (p >= cur_n_) break;
                if (!task.active.test(p)) continue;
                int& degree = task.productive_degree[as_size(p)];
                task.sole_productive_line[as_size(p)] = (degree == 0) ? id : kManyProductiveLines;
                ++degree;
            }
        }
    }

    [[nodiscard]] int find_any_productive_line_for_point(const Task& task, int p) const {
        // ordered_point_lines_[p] is already filtered to activate_at <= cur_n_,
        // so no activate_at guard is needed here.  Lines with activate_at > cur_n_
        // are absent from the list; lines that are blocked or have low coverage
        // have productive_pos == -1 and are skipped by the inner check.
        for (int id : ordered_point_lines_[as_size(p)]) {
            if (task.productive_pos[as_size(id)] >= 0) return id;
        }
        return kNoProductiveLine;
    }

    [[nodiscard]] int exact_unique_points_for_line(const Task& task, int target_line) const {
        int count = 0;
        for (int p : lines_[as_size(target_line)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;

            int exact_sole = kNoProductiveLine;
            int exact_degree = 0;
            for (int id : ordered_point_lines_[as_size(p)]) {
                if (!task.available[as_size(id)]) continue;
                if (task.line_cover[as_size(id)] < 3) continue;
                exact_sole = id;
                if (++exact_degree >= 2) break;
            }
            if (exact_degree == 1 && exact_sole == target_line) ++count;
        }
        return count;
    }

    [[nodiscard]] int exact_unique_points_for_line_with_threshold(const Task& task, int target_line, int threshold) const {
        int clean_unique = 0;
        int dirty_active_points = 0;
        for (int p : lines_[as_size(target_line)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            if (!task.dirty_points[as_size(p)]) {
                if (task.productive_degree[as_size(p)] == 1 && task.sole_productive_line[as_size(p)] == target_line) ++clean_unique;
            } else {
                ++dirty_active_points;
            }
        }
        if (clean_unique + dirty_active_points < threshold) return clean_unique;

        int count = clean_unique;
        int remaining_dirty = dirty_active_points;
        for (int p : lines_[as_size(target_line)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p) || !task.dirty_points[as_size(p)]) continue;

            int exact_sole = kNoProductiveLine;
            int exact_degree = 0;
            for (int id : ordered_point_lines_[as_size(p)]) {
                if (!task.available[as_size(id)] || task.line_cover[as_size(id)] < 3) continue;
                exact_sole = id;
                if (++exact_degree >= 2) break;
            }
            if (exact_degree == 1 && exact_sole == target_line) ++count;
            --remaining_dirty;
            if (count + remaining_dirty < threshold) return count;
        }
        return count;
    }

    [[nodiscard]] ForcedChoice exact_forced_choice(const Task& task, WorkerContext& ctx) const {
        ForcedChoice best;
        auto& exact_unique_count = ctx.forced_unique_count;
        std::vector<int> exact_touched;
        exact_touched.reserve(task.productive_ids.size());

        task.active.for_each_set_bit(words_, [&](int p) {
            int exact_sole = kNoProductiveLine;
            int exact_degree = 0;
            for (int id : ordered_point_lines_[as_size(p)]) {
                // Fast check: is this line currently productive? (available and cover>=3)
                if (task.productive_pos[as_size(id)] < 0) continue;
                exact_sole = id;
                if (++exact_degree >= 2) break;
            }
            if (exact_degree == 1 && exact_unique_count[as_size(exact_sole)]++ == 0)
                exact_touched.push_back(exact_sole);
        });

        for (int id : exact_touched) {
            const int unique_points = exact_unique_count[as_size(id)];
            exact_unique_count[as_size(id)] = 0;
            if (unique_points < 4) continue;

            const int gain = positive_gain(task, id);
            if (gain <= 0) continue;

            const bool better = (unique_points > best.best_unique) ||
                                (unique_points == best.best_unique && gain > best.best_gain) ||
                                (unique_points == best.best_unique && gain == best.best_gain && (best.best_line < 0 || id < best.best_line));
            if (better) { best.best_line = id; best.best_unique = unique_points; best.best_gain = gain; }
        }
        return best;
    }

    void on_productive_line_removed(Task& task, int id) const {
        for (int p : lines_[as_size(id)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;

            int& degree = task.productive_degree[as_size(p)];
            int& sole = task.sole_productive_line[as_size(p)];
            if (degree <= 0) continue;

            if (degree == 1) { degree = 0; sole = kNoProductiveLine; }
            else if (degree == 2) { degree = 1; sole = find_any_productive_line_for_point(task, p); }
            else { --degree; sole = kManyProductiveLines; }
        }
    }

    void on_productive_line_reinserted(Task& task, int id) const {
        for (int p : lines_[as_size(id)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            int& degree = task.productive_degree[as_size(p)];
            if (degree == 0) { degree = 1; task.sole_productive_line[as_size(p)] = id; }
            else { ++degree; task.sole_productive_line[as_size(p)] = kManyProductiveLines; }
        }
    }

    void remove_productive_impl(Task& task, int id, std::vector<int>* removed_log) const {
        int& pos = task.productive_pos[as_size(id)];
        if (pos < 0) return;
        const int last_id = task.productive_ids.back();
        task.productive_ids[as_size(pos)] = last_id;
        task.productive_pos[as_size(last_id)] = pos;
        task.productive_ids.pop_back();
        pos = -1;
        if (removed_log) removed_log->push_back(id);
        on_productive_line_removed(task, id);
    }

    // O(1) swap-to-back removal from productive_ids; records in the worker stack for reversal.
    void remove_productive(Task& task, WorkerContext& ctx, int id) const { remove_productive_impl(task, id, &ctx.removed_from_productive); }

    // Called during undo to re-insert a line that was removed from the productive set.
    void reinsert_productive(Task& task, int id) const {
        if (task.productive_pos[as_size(id)] >= 0) return;
        task.productive_pos[as_size(id)] = as_int(task.productive_ids.size());
        task.productive_ids.push_back(id);
        on_productive_line_reinserted(task, id);
    }

    void add_blocked_productive_impl(Task& task, int id) const {
        int& pos = task.blocked_productive_pos[as_size(id)];
        if (pos >= 0) return;  // already present
        // Insert in sorted order
        auto& vec = task.blocked_productive_ids;
        auto it = std::lower_bound(vec.begin(), vec.end(), id);
        pos = static_cast<int>(it - vec.begin());
        vec.insert(it, id);
        // Update positions of elements shifted right
        for (std::size_t i = pos + 1; i < vec.size(); ++i) {
            task.blocked_productive_pos[as_size(vec[i])] = static_cast<int>(i);
        }
    }

    void remove_blocked_productive_impl(Task& task, int id, std::vector<int>* removed_log) const {
        int& pos = task.blocked_productive_pos[as_size(id)];
        if (pos < 0) return;
        auto& vec = task.blocked_productive_ids;
        // Erase at position pos (the vector is sorted, so we can erase directly)
        vec.erase(vec.begin() + pos);
        // Update positions of elements after the erased one
        for (std::size_t i = pos; i < vec.size(); ++i) {
            task.blocked_productive_pos[as_size(vec[i])] = static_cast<int>(i);
        }
        pos = -1;
        if (removed_log) removed_log->push_back(id);
    }

    void remove_blocked_productive(Task& task, WorkerContext& ctx, int id) const {
        remove_blocked_productive_impl(task, id, &ctx.removed_from_blocked_productive);
    }

    void disable_line_without_undo(Task& task, int id) const {
        if (!task.available[as_size(id)]) return;
        task.available[as_size(id)] = 0;
        if (task.line_cover[as_size(id)] >= 3) add_blocked_productive_impl(task, id);
        remove_productive_impl(task, id, nullptr);
    }

    void disable_line(Task& task, WorkerContext& ctx, int id) const {
        if (!task.available[as_size(id)]) return;
        task.available[as_size(id)] = 0;
        if (task.line_cover[as_size(id)] >= 3) add_blocked_productive_impl(task, id);
        ctx.disabled_lines.push_back(id);
        remove_productive(task, ctx, id);  // unavailable => non-productive
    }

    void cover_point(Task& task, WorkerContext& ctx, int p) const {
        if (!task.active.test(p)) return;
        ctx.covered_point_degree.push_back(task.productive_degree[as_size(p)]);
        ctx.covered_point_sole.push_back(task.sole_productive_line[as_size(p)]);
        ctx.covered_point_dirty.push_back(task.dirty_points[as_size(p)]);
        task.productive_degree[as_size(p)] = 0;
        task.sole_productive_line[as_size(p)] = kNoProductiveLine;
        task.dirty_points[as_size(p)] = 0;
        task.active.reset(p);
        --task.active_count;
        ctx.newly_covered_points.push_back(p);
        for (int id : ordered_point_lines_[as_size(p)]) {
            const bool available = task.available[as_size(id)] != 0;
            const bool blocked_productive = task.blocked_productive_pos[as_size(id)] >= 0;
            if (!available && !blocked_productive) continue;
            const uint16_t old_cover = task.line_cover[as_size(id)]--;
            ctx.decremented_lines.push_back(id);
            if (old_cover != 3) continue;
            if (available) remove_productive(task, ctx, id);  // cover 3->2: no longer productive
            else remove_blocked_productive(task, ctx, id);
        }
    }

    void select_line(Task& task, WorkerContext& ctx, int id) const {
        disable_line(task, ctx, id);
        for (int p : lines_[as_size(id)].points) {
            if (p >= cur_n_) break;
            cover_point(task, ctx, p);
        }
    }

    void undo_to(Task& task, WorkerContext& ctx, UndoMark mark) const {
        for (std::size_t i = mark.decremented_lines; i < ctx.decremented_lines.size(); ++i) ++task.line_cover[as_size(ctx.decremented_lines[i])];
        for (std::size_t i = mark.newly_covered_points; i < ctx.newly_covered_points.size(); ++i) {
            task.active.set(ctx.newly_covered_points[i]);
            ++task.active_count;
        }
        for (std::size_t i = mark.removed_from_blocked_productive; i < ctx.removed_from_blocked_productive.size(); ++i)
            add_blocked_productive_impl(task, ctx.removed_from_blocked_productive[i]);
        for (std::size_t i = mark.disabled_lines; i < ctx.disabled_lines.size(); ++i) {
            remove_blocked_productive_impl(task, ctx.disabled_lines[i], nullptr);
            task.available[as_size(ctx.disabled_lines[i])] = 1;
        }
        for (std::size_t i = mark.removed_from_productive; i < ctx.removed_from_productive.size(); ++i) reinsert_productive(task, ctx.removed_from_productive[i]);
        for (std::size_t i = mark.covered_point_metadata; i < ctx.covered_point_degree.size(); ++i) {
            const int p = ctx.newly_covered_points[i];
            task.productive_degree[as_size(p)] = ctx.covered_point_degree[i];
            task.sole_productive_line[as_size(p)] = ctx.covered_point_sole[i];
            task.dirty_points[as_size(p)] = ctx.covered_point_dirty[i];
        }
        ctx.removed_from_productive.resize(mark.removed_from_productive);
        ctx.removed_from_blocked_productive.resize(mark.removed_from_blocked_productive);
        ctx.decremented_lines.resize(mark.decremented_lines);
        ctx.newly_covered_points.resize(mark.newly_covered_points);
        ctx.covered_point_degree.resize(mark.covered_point_metadata);
        ctx.covered_point_sole.resize(mark.covered_point_metadata);
        ctx.covered_point_dirty.resize(mark.covered_point_metadata);
        ctx.disabled_lines.resize(mark.disabled_lines);
    }

    void select_line_without_undo(Task& task, int id) const {
        if (!task.available[as_size(id)]) return;
        disable_line_without_undo(task, id);

        for (int p : lines_[as_size(id)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            task.dirty_points[as_size(p)] = 0;
            task.active.reset(p);
            --task.active_count;
            for (int lid : ordered_point_lines_[as_size(p)]) {
                const bool available = task.available[as_size(lid)] != 0;
                const bool blocked_productive = task.blocked_productive_pos[as_size(lid)] >= 0;
                if (!available && !blocked_productive) continue;
                if (task.line_cover[as_size(lid)]-- != 3) continue;
                if (available) remove_productive_impl(task, lid, nullptr);  // cover 3->2: drop from productive set
                else remove_blocked_productive_impl(task, lid, nullptr);
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
            if (id < 0 || id >= as_int(lines_.size()) || lines_[as_size(id)].activate_at > cur_n_) continue;
            const int delta = positive_gain(trial, id);
            if (delta <= 0) continue;
            gain += delta;
            chosen.push_back(id);
            select_line_without_undo(trial, id);
        }
        if (gain > current_best_gain()) submit_candidate(gain, chosen);
    }

    void greedy_complete_from_seed(const Task& base, WorkerContext& ctx, const std::vector<int>& line_ids) {
        if (line_ids.empty()) return;
        Task trial = base;
        ensure_task_search_buffers(trial);
        for (int id : line_ids) {
            if (id < 0 || id >= as_int(lines_.size()) || lines_[as_size(id)].activate_at > cur_n_) continue;
            const int delta = positive_gain(trial, id);
            if (delta <= 0) continue;
            select_line_without_undo(trial, id);
            trial.current_gain += delta;
            trial.current_choice.push_back(id);
        }
        maybe_submit_greedy(trial, ctx);
    }

    [[nodiscard]] static bool better_branch_candidate(
        int gain_a, long double slack_a, long long pressure_a, int id_a,
        int gain_b, long double slack_b, long long pressure_b, int id_b) noexcept
    {
        return (gain_a > gain_b) ||
               (gain_a == gain_b && slack_a > slack_b + 1e-18L) ||
               (gain_a == gain_b && std::abs(slack_a - slack_b) <= 1e-18L && pressure_a > pressure_b) ||
               (gain_a == gain_b && std::abs(slack_a - slack_b) <= 1e-18L && pressure_a == pressure_b && id_a < id_b);
    }

    [[nodiscard]] BranchChoice choose_branch(Task& task, WorkerContext& ctx, long double node_gain_ub) const {
        struct Candidate {
            int         id = -1;
            int         gain = -1;
            long double slack = -1.0L;
            long long   pressure = -1;
            long double exclude_ub = 0.0L;
        };

        std::array<Candidate, kStrongBranchShortlist> shortlist{};
        int shortlist_size = 0;

        auto branch_pressure = [&](Candidate& candidate) -> long long {
            if (candidate.pressure >= 0) return candidate.pressure;
            long long pressure = 0;
            for (int p : lines_[as_size(candidate.id)].points) {
                if (p >= cur_n_) break;
                if (task.active.test(p)) pressure += static_cast<long long>(ordered_point_lines_[as_size(p)].size());
            }
            return candidate.pressure = pressure;
        };

        auto better_shortlist_candidate = [&](Candidate& candidate_a, Candidate& candidate_b) -> bool {
            if (candidate_a.gain != candidate_b.gain) return candidate_a.gain > candidate_b.gain;
            if (candidate_a.slack > candidate_b.slack + 1e-18L) return true;
            if (candidate_b.slack > candidate_a.slack + 1e-18L) return false;
            const long long pressure_a = branch_pressure(candidate_a);
            const long long pressure_b = branch_pressure(candidate_b);
            return (pressure_a > pressure_b) || (pressure_a == pressure_b && candidate_a.id < candidate_b.id);
        };

        auto insert_candidate = [&](Candidate candidate) {
            int pos = shortlist_size;
            if (pos < kStrongBranchShortlist) { shortlist[as_size(pos)] = candidate; ++shortlist_size; } 
            else if (!better_shortlist_candidate(candidate, shortlist.back())) return;
            else { shortlist.back() = candidate; pos = kStrongBranchShortlist - 1; }

            while (pos > 0 && better_shortlist_candidate(shortlist[as_size(pos)], shortlist[as_size(pos - 1)])) {
                std::swap(shortlist[as_size(pos)], shortlist[as_size(pos - 1)]);
                --pos;
            }
        };

        const auto& line_slack = ctx.lagrangian.line_slack;
        for (int id : task.productive_ids) {
            insert_candidate({ id, as_int(task.line_cover[as_size(id)]) - 2, std::max(0.0L, line_slack[as_size(id)]), -1, 0.0L });
        }

        if (shortlist_size == 0) return {};

        for (int i = 0; i < shortlist_size; ++i) {
            shortlist[as_size(i)].exclude_ub = exclude_branch_upper_bound(task, ctx, shortlist[as_size(i)].id, node_gain_ub);
        }

        BranchChoice best{ shortlist[0].id, shortlist[0].exclude_ub };

        const int node_cost_lb = cost_from_gain(cur_n_, as_int(std::floor(node_gain_ub + 1e-9L)));
        const bool use_strong_branching = shortlist_size > 1 && current_best_cost() - node_cost_lb <= 1 &&
                                          task.active_count >= 64 && task.current_choice.size() <= kStrongBranchMaxDepth;
        if (!use_strong_branching) return best;
        strong_branch_atomic_ptr->fetch_add(1, std::memory_order_relaxed);

        const std::vector<int> parent_active_points = ctx.lagrangian.active_points;
        int best_pruned_children = -1;
        long double best_worst_child_ub = std::numeric_limits<long double>::infinity();
        long double best_total_child_ub = std::numeric_limits<long double>::infinity();
        int best_gain = -1;
        long double best_slack = -1.0L;
        long long best_pressure = -1;

        for (int i = 0; i < shortlist_size; ++i) {
            Candidate& candidate = shortlist[as_size(i)];
            const int start_gain = task.current_gain;
            const std::size_t start_choice_size = task.current_choice.size();
            const UndoMark undo = mark_undo(ctx);
            const std::size_t dual_mark = save_dual_for_points(ctx, parent_active_points);

            long double include_ub = static_cast<long double>(task.current_gain);
            const int delta = positive_gain(task, candidate.id);
            if (delta > 0) {
                select_line(task, ctx, candidate.id);
                task.current_gain += delta;
                task.current_choice.push_back(candidate.id);
                apply_forced_line_reduction(task, ctx, true);
                include_ub = lagrangian_upper_bound(task, ctx, kStrongBranchMaxIters, nullptr);
            }

            restore_saved_dual(ctx, dual_mark);
            task.current_choice.resize(start_choice_size);
            task.current_gain = start_gain;
            undo_to(task, ctx, undo);

            const int pruned_children = as_int(gain_upper_bound_cannot_beat_best_cost(include_ub)) +
                                        as_int(gain_upper_bound_cannot_beat_best_cost(candidate.exclude_ub));
            const long double worst_child_ub = std::max(include_ub, candidate.exclude_ub);
            const long double total_child_ub = include_ub + candidate.exclude_ub;
            const long long candidate_pressure = branch_pressure(candidate);

            const bool better = (pruned_children > best_pruned_children) ||
                                (pruned_children == best_pruned_children && worst_child_ub + 1e-9L < best_worst_child_ub) ||
                                (pruned_children == best_pruned_children && std::abs(worst_child_ub - best_worst_child_ub) <= 1e-9L && total_child_ub + 1e-9L < best_total_child_ub) ||
                                (pruned_children == best_pruned_children && std::abs(worst_child_ub - best_worst_child_ub) <= 1e-9L && std::abs(total_child_ub - best_total_child_ub) <= 1e-9L &&
                                 better_branch_candidate(candidate.gain, candidate.slack, candidate_pressure, candidate.id, best_gain, best_slack, best_pressure, best.line));
            if (better) {
                best_pruned_children = pruned_children; best_worst_child_ub = worst_child_ub; best_total_child_ub = total_child_ub;
                best_gain = candidate.gain; best_slack = candidate.slack; best_pressure = candidate_pressure;
                best.line = candidate.id; best.exclude_ub = candidate.exclude_ub;
            }
        }
        return best;
    }

    [[nodiscard]] static std::pair<int, int> root_productive_metrics(const Task& task) noexcept {
        int productive_incidence = 0;
        for (int id : task.productive_ids) productive_incidence += as_int(task.line_cover[as_size(id)]);
        return { as_int(task.productive_ids.size()), productive_incidence };
    }

    void apply_forced_line_reduction(Task& task, WorkerContext& ctx, bool record_undo) const {
        auto& unique_count = ctx.forced_unique_count;
        auto& dirty_unique_count = ctx.forced_dirty_unique_count;
        auto& touched = ctx.forced_touched_lines;

        while (true) {
            touched.clear();
            task.active.for_each_set_bit(words_, [&](int p) {
                const int unique_line = task.sole_productive_line[as_size(p)];
                if (unique_line >= 0) {
                    if (unique_count[as_size(unique_line)]++ == 0) touched.push_back(unique_line);
                    dirty_unique_count[as_size(unique_line)] += as_int(task.dirty_points[as_size(p)] != 0);
                }
            });

            int best_line = -1, best_unique = 2, best_gain = -1;
            for (int id : touched) {
                const int unique_points = unique_count[as_size(id)];
                const int dirty_unique_points = dirty_unique_count[as_size(id)];
                unique_count[as_size(id)] = 0;
                dirty_unique_count[as_size(id)] = 0;
                if (unique_points < 4) continue;

                const int gain = positive_gain(task, id);
                if (gain <= 0) continue;

                int checked_unique = unique_points;
                if (dirty_unique_points > 0 && unique_points - dirty_unique_points < 4) {
                    checked_unique = exact_unique_points_for_line_with_threshold(task, id, 4);
                    if (checked_unique < 4) continue;
                }

                const bool better = (checked_unique > best_unique) || (checked_unique == best_unique && gain > best_gain) ||
                                    (checked_unique == best_unique && gain == best_gain && (best_line < 0 || id < best_line));
                if (better) { best_line = id; best_unique = checked_unique; best_gain = gain; }
            }

            if (best_line < 0) break;

            if (record_undo) {
                const int exact_unique = exact_unique_points_for_line(task, best_line);
                if (exact_unique < 4) {
                    // Metadata overestimated; abort forced reduction at this node.
                    break;
                }
                best_unique = exact_unique;
                best_gain = positive_gain(task, best_line);
            }

            if (trace_this_n_ && record_undo) {
                const int exact_meta_unique = exact_unique_points_for_line(task, best_line);
                const ForcedChoice exact = exact_forced_choice(task, ctx);
                if (exact_meta_unique != best_unique || exact.best_line != best_line) {
                    debug_append_text(".events.txt", std::format("FALSE_FORCE depth={} chosen={} best_line={} meta_unique={} meta_gain={} exact_meta_unique={} exact_line={} exact_unique={} exact_gain={}\n",
                        task.current_choice.size(), task.current_choice.size(), best_line, best_unique, best_gain, exact_meta_unique, exact.best_line, exact.best_unique, exact.best_gain));
                }
                if (force_exact_forced_lines_) { best_line = exact.best_line; best_unique = exact.best_unique; best_gain = exact.best_gain; }
            }

            if (best_line < 0) break;

            // Safe kernel rule: if a productive heavy line owns at least four
            // currently active points that no other productive heavy line can
            // cover, there is always an optimal completion that includes it.
            // Note: the threshold of four (not three) is intentional — lowering
            // it to three was tested over N=730..800 and produced incorrect
            // optimal line counts at N=731, confirming it is not a safe forcing rule.
            forced_atomic_ptr->fetch_add(1, std::memory_order_relaxed);
            if (record_undo) select_line(task, ctx, best_line);
            else select_line_without_undo(task, best_line);
            task.current_gain += best_gain;
            task.current_choice.push_back(best_line);
        }
    }

    [[nodiscard]] int greedy_completion(const Task& task, WorkerContext& ctx, bool keep_lines) const {
        BitMask1024 temp_active = task.active;
        const int max_gain = std::max(0, task.active_count - 2);
        std::fill_n(ctx.greedy_bucket_head.begin(), max_gain + 1, -1);
        for (int id : task.productive_ids) {
            const int gain = as_int(task.line_cover[as_size(id)]) - 2;
            if (gain <= 0) continue;
            ctx.greedy_bucket_next[as_size(id)] = ctx.greedy_bucket_head[as_size(gain)];
            ctx.greedy_bucket_head[as_size(gain)] = id;
        }

        auto& chosen = ctx.greedy_lines;
        chosen.clear();

        int gain = 0;
        for (int bucket = max_gain; bucket > 0; --bucket) {
            for (int id = ctx.greedy_bucket_head[as_size(bucket)]; id >= 0; id = ctx.greedy_bucket_next[as_size(id)]) {
                int current_gain = -2;
                for (int p : lines_[as_size(id)].points) {
                    if (p >= cur_n_) break;
                    if (temp_active.test(p)) ++current_gain;
                }
                if (current_gain <= 0) continue;

                gain += current_gain;
                if (keep_lines) chosen.push_back(id);
                for (int p : lines_[as_size(id)].points) {
                    if (p >= cur_n_) break;
                    if (temp_active.test(p)) temp_active.reset(p);
                }
            }
        }
        return gain;
    }

    [[nodiscard]] int greedy_completion_compact(
        const Task& task, WorkerContext& ctx, const std::vector<int>& compact_line_ids,
        const std::vector<std::size_t>& compact_offsets, const std::vector<CompactPointIndex>& compact_points, bool keep_lines) const
    {
        BitMask1024 temp_active = task.active;
        const int max_gain = std::max(0, task.active_count - 2);
        std::fill_n(ctx.greedy_bucket_head.begin(), max_gain + 1, -1);
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            const int line_gain = as_int(task.line_cover[as_size(compact_line_ids[li])]) - 2;
            if (line_gain <= 0) continue;
            ctx.greedy_bucket_next[li] = ctx.greedy_bucket_head[as_size(line_gain)];
            ctx.greedy_bucket_head[as_size(line_gain)] = as_int(li);
        }

        auto& chosen = ctx.greedy_lines;
        chosen.clear();

        int gain = 0;
        const auto* points = compact_points.data();
        for (int bucket = max_gain; bucket > 0; --bucket) {
            for (int li = ctx.greedy_bucket_head[as_size(bucket)]; li >= 0; li = ctx.greedy_bucket_next[as_size(li)]) {
                const std::size_t line_index = as_size(li);
                int current_gain = -2;
                for (std::size_t pos = compact_offsets[line_index]; pos < compact_offsets[line_index + 1]; ++pos) {
                    if (temp_active.test(static_cast<int>(points[pos]))) ++current_gain;
                }
                if (current_gain <= 0) continue;

                gain += current_gain;
                if (keep_lines) chosen.push_back(compact_line_ids[line_index]);
                for (std::size_t pos = compact_offsets[line_index]; pos < compact_offsets[line_index + 1]; ++pos) {
                    const int point = static_cast<int>(points[pos]);
                    if (temp_active.test(point)) temp_active.reset(point);
                }
            }
        }
        return gain;
    }

    void build_compact_productive_lines(const Task& task, LagrangianScratch& scratch) const {
        auto& active_points = scratch.active_points;
        auto& compact_line_ids = scratch.compact_line_ids;
        auto& compact_offsets = scratch.compact_line_offsets;
        auto& compact_points = scratch.compact_line_points;
        auto& compact_point_aps = scratch.compact_point_aps;
        auto& compact_eval3 = scratch.compact_eval3;
        auto& compact_eval4 = scratch.compact_eval4;
        auto& compact_eval5 = scratch.compact_eval5;
        auto& compact_eval_other_lines = scratch.compact_eval_other_lines;

        active_points.clear();
        compact_line_ids.clear(); compact_offsets.clear(); compact_points.clear(); compact_point_aps.clear();
        compact_eval3.clear(); compact_eval4.clear(); compact_eval5.clear(); compact_eval_other_lines.clear();
        const std::size_t productive_count = task.productive_ids.size();
        
        if (compact_line_ids.capacity() < productive_count) compact_line_ids.reserve(productive_count);
        if (compact_offsets.capacity() < productive_count + 1) compact_offsets.reserve(productive_count + 1);

        std::size_t total_active_incidence = 0;
        for (int id : task.productive_ids) total_active_incidence += as_size(task.line_cover[as_size(id)]);
        if (compact_points.capacity() < total_active_incidence) compact_points.reserve(total_active_incidence);
        if (compact_point_aps.capacity() < total_active_incidence) compact_point_aps.reserve(total_active_incidence);

        BitMask1024 productive_union;
        compact_offsets.push_back(0);
        for (int id : task.productive_ids) {
            compact_line_ids.push_back(id);
            int remaining = as_int(task.line_cover[as_size(id)]);
            for (int p : lines_[as_size(id)].points) {
                if (p >= cur_n_) break;
                if (!task.active.test(p)) continue;
                productive_union.set(p);
                compact_points.push_back(static_cast<CompactPointIndex>(p));
                if (--remaining == 0) break;
            }
            compact_offsets.push_back(compact_points.size());
        }

        // Reuse the compact productive-line pass to recover the active-point set
        // for this node, instead of rescanning every active point through the
        // global incidence lists first and rebuilding the same support again here.
        productive_union.for_each_set_bit(words_, [&](int p) { active_points.push_back(p); });

        {
            // Precompute (line_cover[id] - 2) for every compact line.
            // task.line_cover is ~16 KB and is accessed with scattered IDs inside
            // eval_dense_fast / eval_dense_exact / the polish initialiser; each of
            // those hot paths is called up to 64+ times per bound invocation.
            // Replacing those random loads with sequential reads from this small
            // compact array reduces the L1 working set from ~34 KB to ~18 KB,
            // eliminating cache-eviction pressure on the other hot arrays (yy,
            // compact_points, subgrad).
            auto& rem2 = scratch.compact_rem2;
            rem2.resize(compact_line_ids.size());
            for (std::size_t li = 0; li < compact_line_ids.size(); ++li)
                rem2[li] = as_int(task.line_cover[as_size(compact_line_ids[li])]) - 2;
        }

        // -----------------------------------------------------------------
        // Build compact inverse: ap index → productive-line compact indices.
        //
        // The coordinate polish iterates active_points and, for each point p,
        // needs all productive lines covering p.  Scanning incidence_[p] and
        // filtering is wasteful at deeper DFS nodes where many lines are
        // unavailable or have coverage < 3.  The compact inverse contains only
        // productive entries (guaranteed by the productive invariant), so no
        // filtering is needed and the total scan is bounded by pinc rather
        // than global incidence.
        //
        // Build cost: O(ap_size + pinc).  Reset cost: O(ap_size).
        // -----------------------------------------------------------------
        const int ap_size = as_int(active_points.size());
        auto& ap_offsets  = scratch.compact_pt_prod_offsets;
        auto& ap_flat     = scratch.compact_pt_prod_flat;
        auto& fill_pos    = scratch.compact_ap_fill_pos;
        auto& pt_to_ap    = scratch.point_to_ap_scratch;

        // Map global point index → position in active_points.
        // pt_to_ap is kept all-(-1) outside this function (invariant).
        for (int ai = 0; ai < ap_size; ++ai) pt_to_ap[as_size(active_points[ai])] = ai;

        // Count pass: tally productive lines per ap index (stored at ap_offsets[ai+1]).
        ap_offsets.assign(as_size(ap_size + 1), 0);
        compact_point_aps.resize(compact_points.size());
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            for (std::size_t pos = compact_offsets[li]; pos < compact_offsets[li + 1]; ++pos) {
                const int ai = pt_to_ap[as_size(compact_points[pos])];
                compact_point_aps[pos] = static_cast<CompactPointIndex>(ai);
                ++ap_offsets[as_size(ai + 1)];
            }
        }

        compact_eval3.reserve(compact_line_ids.size());
        compact_eval4.reserve(compact_line_ids.size());
        compact_eval5.reserve(compact_line_ids.size());
        compact_eval_other_lines.reserve(compact_line_ids.size());
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            const std::size_t begin = compact_offsets[li];
            const std::size_t end = compact_offsets[li + 1];
            const std::size_t rem = end - begin;
            if (rem == 3) {
                compact_eval3.push_back({
                    compact_point_aps[begin],
                    compact_point_aps[begin + 1],
                    compact_point_aps[begin + 2]
                });
            } else if (rem == 4) {
                compact_eval4.push_back({
                    compact_point_aps[begin],
                    compact_point_aps[begin + 1],
                    compact_point_aps[begin + 2],
                    compact_point_aps[begin + 3]
                });
            } else if (rem == 5) {
                compact_eval5.push_back({
                    compact_point_aps[begin],
                    compact_point_aps[begin + 1],
                    compact_point_aps[begin + 2],
                    compact_point_aps[begin + 3],
                    compact_point_aps[begin + 4]
                });
            } else {
                compact_eval_other_lines.push_back(as_int(li));
            }
        }

        // Prefix-sum to get start offsets.
        for (int i = 0; i < ap_size; ++i) ap_offsets[as_size(i + 1)] += ap_offsets[as_size(i)];

        // Allocate flat buffer (total entries = pinc).
        ap_flat.resize(as_size(ap_offsets[ap_size]));

        // Copy start offsets as write cursors.
        fill_pos.resize(as_size(ap_size));
        for (int i = 0; i < ap_size; ++i) fill_pos[as_size(i)] = ap_offsets[as_size(i)];

        // Fill pass: write compact line index li for every (ap, li) pair.
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            for (std::size_t pos = compact_offsets[li]; pos < compact_offsets[li + 1]; ++pos) {
                const int ai = compact_point_aps[pos];
                ap_flat[as_size(fill_pos[as_size(ai)]++)] = as_int(li);
            }
        }

        // Restore pt_to_ap invariant: reset only the entries we touched.
        for (int p : active_points) pt_to_ap[as_size(p)] = -1;

        // Size the compact slack array (will be filled in the polish setup).
        scratch.work_line_slack_compact.resize(compact_line_ids.size());
    }

    [[nodiscard]] long double lagrangian_upper_bound(const Task& task, WorkerContext& ctx, int max_iters, int* greedy_gain_out) const {
        // Zero-overhead cancellation: one relaxed atomic load, amortised over all the
        // expensive work below (greedy_completion + 36 subgradient iters).  Returning
        // -1 makes gain_upper_bound_cannot_beat_best_cost() prune this node immediately
        // in dfs() without touching any other hot-path code.
        if (abort_.load(std::memory_order_relaxed)) return -1.0L;

        auto& scratch = ctx.lagrangian;
        auto& active_y = scratch.active_y;
        auto& active_best_y = scratch.active_best_y;
        auto& best_y = scratch.best_y;
        auto& active_subgrad = scratch.active_subgrad;
        auto& line_slack = scratch.line_slack;

        build_compact_productive_lines(task, scratch);
        const auto& active_points = scratch.active_points;
        const auto& compact_line_ids = scratch.compact_line_ids;
        const auto& compact_offsets = scratch.compact_line_offsets;
        const auto& compact_points = scratch.compact_line_points;
        const auto& compact_point_aps = scratch.compact_point_aps;
        const auto& compact_eval3 = scratch.compact_eval3;
        const auto& compact_eval4 = scratch.compact_eval4;
        const auto& compact_eval5 = scratch.compact_eval5;
        const auto& compact_eval_other_lines = scratch.compact_eval_other_lines;
        const auto& compact_rem2 = scratch.compact_rem2;
        const auto& ap_offsets = scratch.compact_pt_prod_offsets;
        const auto& ap_flat = scratch.compact_pt_prod_flat;
        const int active_point_count = as_int(active_points.size());
        const auto* compact_offsets_data = compact_offsets.data();
        const auto* compact_point_aps_data = compact_point_aps.data();
        const auto* compact_rem2_data = compact_rem2.data();

        if (active_points.empty()) {
            scratch.dual_ready = false;
            return static_cast<long double>(task.current_gain);
        }

        if (scratch.dual_ready) {
            for (int ai = 0; ai < active_point_count; ++ai) {
                const int p = active_points[as_size(ai)];
                active_y[as_size(ai)] = std::clamp<DualValue>(best_y[as_size(p)], 0.0, 1.0);
            }
        } else {
            for (int ai = 0; ai < active_point_count; ++ai) {
                DualValue best_density = 0.0;
                for (int k = ap_offsets[as_size(ai)]; k < ap_offsets[as_size(ai + 1)]; ++k) {
                    const int li = ap_flat[as_size(k)];
                    const int rem2 = compact_rem2[as_size(li)];
                    const int rem = rem2 + 2;
                    best_density = std::max(best_density, static_cast<DualValue>(rem2) / static_cast<DualValue>(rem));
                }
                active_y[as_size(ai)] = std::min<DualValue>(1.0, best_density);
            }
        }

        auto eval_dense_fast = [&](const std::vector<DualValue>& yy, std::vector<DualValue>* subgradient) -> DualValue {
            const auto* yy_data = yy.data();
            auto* subgradient_data = subgradient ? subgradient->data() : nullptr;
            DualValue value = 0.0;
            for (int ai = 0; ai < active_point_count; ++ai) {
                value += yy_data[ai];
                if (subgradient_data) subgradient_data[ai] = 1.0;
            }
            for (const CompactLine3& line : compact_eval3) {
                const DualValue rc =
                    1.0 - (yy_data[line.a] + yy_data[line.b] + yy_data[line.c]);
                if (rc <= 1e-18) continue;
                value += rc;
                if (subgradient_data) {
                    subgradient_data[line.a] -= 1.0;
                    subgradient_data[line.b] -= 1.0;
                    subgradient_data[line.c] -= 1.0;
                }
            }
            for (const CompactLine4& line : compact_eval4) {
                const DualValue rc =
                    2.0 - (yy_data[line.a] + yy_data[line.b] + yy_data[line.c] + yy_data[line.d]);
                if (rc <= 1e-18) continue;
                value += rc;
                if (subgradient_data) {
                    subgradient_data[line.a] -= 1.0;
                    subgradient_data[line.b] -= 1.0;
                    subgradient_data[line.c] -= 1.0;
                    subgradient_data[line.d] -= 1.0;
                }
            }
            for (const CompactLine5& line : compact_eval5) {
                const DualValue rc =
                    3.0 - (yy_data[line.a] + yy_data[line.b] + yy_data[line.c] + yy_data[line.d] + yy_data[line.e]);
                if (rc <= 1e-18) continue;
                value += rc;
                if (subgradient_data) {
                    subgradient_data[line.a] -= 1.0;
                    subgradient_data[line.b] -= 1.0;
                    subgradient_data[line.c] -= 1.0;
                    subgradient_data[line.d] -= 1.0;
                    subgradient_data[line.e] -= 1.0;
                }
            }
            for (int li_int : compact_eval_other_lines) {
                const std::size_t li = as_size(li_int);
                const std::size_t begin = compact_offsets_data[li];
                const std::size_t end = compact_offsets_data[li + 1];
                DualValue sum_y = 0.0;
                for (std::size_t pos = begin; pos < end; ++pos)
                    sum_y += yy_data[compact_point_aps_data[pos]];
                const DualValue rc = static_cast<DualValue>(compact_rem2_data[li]) - sum_y;
                if (rc > 1e-18) {
                    value += rc;
                    if (subgradient_data) {
                        for (std::size_t pos = begin; pos < end; ++pos)
                            subgradient_data[compact_point_aps_data[pos]] -= 1.0;
                    }
                }
            }
            return value;
        };

        auto eval_dense_exact = [&](const std::vector<DualValue>& yy, bool store_slack) -> long double {
            const auto* yy_data = yy.data();
            long double value = 0.0L;
            for (int ai = 0; ai < active_point_count; ++ai)
                value += static_cast<long double>(yy_data[ai]);
            for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
                const int id = compact_line_ids[li];
                const std::size_t begin = compact_offsets_data[li];
                const std::size_t end = compact_offsets_data[li + 1];
                long double sum_y = 0.0L;
                for (std::size_t pos = begin; pos < end; ++pos)
                    sum_y += static_cast<long double>(yy_data[compact_point_aps_data[pos]]);
                const long double rc = static_cast<long double>(compact_rem2_data[li]) - sum_y;
                if (store_slack) line_slack[as_size(id)] = rc;
                if (rc > 1e-18L) value += rc;
            }
            return value;
        };

        auto store_best_dual = [&]() {
            for (int ai = 0; ai < active_point_count; ++ai)
                best_y[as_size(active_points[as_size(ai)])] = active_best_y[as_size(ai)];
        };

        const long double best_exact = eval_dense_exact(active_y, false);
        DualValue best = static_cast<DualValue>(best_exact);
        for (int ai = 0; ai < active_point_count; ++ai)
            active_best_y[as_size(ai)] = active_y[as_size(ai)];

        const long double initial_total = static_cast<long double>(task.current_gain) + best_exact;
        scratch.dual_ready = true;
        if (gain_upper_bound_cannot_beat_best_cost(initial_total)) {
            store_best_dual();
            return initial_total;
        }

        const bool keep_greedy_lines = (greedy_gain_out != nullptr);
        const int greedy_lb = greedy_completion_compact(task, ctx, compact_line_ids, compact_offsets, compact_points, keep_greedy_lines);
        if (greedy_gain_out) *greedy_gain_out = greedy_lb;

        long double confirmed_prune_total = 0.0L;
        auto maybe_confirm_prune = [&](DualValue candidate_best) -> bool {
            const long double optimistic_total = static_cast<long double>(task.current_gain) + static_cast<long double>(candidate_best);
            if (optimistic_total + 1e-9L >= static_cast<long double>(gain_needed_to_beat_best_cost())) return false;
            const long double exact_total = static_cast<long double>(task.current_gain) + eval_dense_exact(active_best_y, false);
            if (!gain_upper_bound_cannot_beat_best_cost(exact_total)) return false;
            confirmed_prune_total = exact_total;
            store_best_dual();
            return true;
        };

        for (int iter = 0; iter < max_iters; ++iter) {
            lag_iters_atomic_ptr->fetch_add(1, std::memory_order_relaxed);
            const DualValue value = eval_dense_fast(active_y, &active_subgrad);
            if (value + 1e-18 < best) {
                best = value;
                for (int ai = 0; ai < active_point_count; ++ai)
                    active_best_y[as_size(ai)] = active_y[as_size(ai)];
                if (maybe_confirm_prune(best)) return confirmed_prune_total;
            }

            DualValue norm2 = 0.0;
            for (int ai = 0; ai < active_point_count; ++ai) {
                const DualValue g = active_subgrad[as_size(ai)];
                norm2 += g * g;
            }
            if (norm2 <= 1e-30) break;

            const DualValue incumbent_gap = static_cast<DualValue>(std::max(0, gain_needed_to_beat_best_cost() - task.current_gain));
            const DualValue target = std::max<DualValue>(static_cast<DualValue>(greedy_lb), incumbent_gap);
            const DualValue gap = std::max<DualValue>(0.0, value - target);
            if (gap <= 1e-15) break;

            const DualValue step = (1.35 * gap) / norm2;
            for (int ai = 0; ai < active_point_count; ++ai) {
                DualValue next = active_y[as_size(ai)] - step * active_subgrad[as_size(ai)];
                if (next < 0.0) next = 0.0;
                else if (next > 1.0) next = 1.0;
                active_y[as_size(ai)] = next;
            }
        }

        if (maybe_confirm_prune(best)) return confirmed_prune_total;

        constexpr int kCoordinatePolishMaxSweeps = 6;
        constexpr long double kCoordinatePolishWindow = 2.0L;
        const long double best_total = static_cast<long double>(task.current_gain) + static_cast<long double>(best);
        if (kCoordinatePolishMaxSweeps > 0 &&
            best_total < static_cast<long double>(gain_needed_to_beat_best_cost()) + kCoordinatePolishWindow) {
            // Initialise compact slack: work_line_slack_compact[li] = (rem-2) - sum_y over
            // active points of line li.  Indexed by compact line index for sequential access.
            auto& wls_c = scratch.work_line_slack_compact;
            for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
                const std::size_t begin = compact_offsets_data[li];
                const std::size_t end = compact_offsets_data[li + 1];
                DualValue slack = static_cast<DualValue>(compact_rem2_data[li]);
                for (std::size_t pos = begin; pos < end; ++pos)
                    slack -= active_best_y[compact_point_aps_data[pos]];
                wls_c[li] = slack;
            }

            auto slack_objective = [&]() -> DualValue {
                DualValue value = 0.0;
                for (int ai = 0; ai < active_point_count; ++ai)
                    value += active_best_y[as_size(ai)];
                for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
                    if (wls_c[li] > 1e-18) value += wls_c[li];
                }
                return value;
            };

            DualValue polished_best = best;
            for (int sweep = 0; sweep < kCoordinatePolishMaxSweeps; ++sweep) {
                bool changed = false;
                // Iterate by ap index so we can address ap_offsets[ai] directly.
                for (int ai = 0; ai < active_point_count; ++ai) {
                    const DualValue old = active_best_y[as_size(ai)];
                    DualValue largest = 0.0;
                    DualValue second  = 0.0;
                    const int beg = ap_offsets[as_size(ai)];
                    const int end = ap_offsets[as_size(ai + 1)];
                    // Scan only productive lines covering p — no availability or coverage
                    // checks needed because the compact invariant guarantees all entries
                    // satisfy available && line_cover >= 3.
                    for (int k = beg; k < end; ++k) {
                        const int li = ap_flat[as_size(k)];
                        const DualValue adjusted = wls_c[as_size(li)] + old;
                        if (adjusted <= 1e-18) continue;
                        if (adjusted >= largest) { second = largest; largest = adjusted; } 
                        else if (adjusted > second) { second = adjusted; }
                    }

                    const DualValue next = std::min<DualValue>(1.0, second);
                    if (std::abs(next - old) <= 1e-15) continue;
                    changed = true;
                    active_best_y[as_size(ai)] = next;
                    const DualValue delta = next - old;
                    for (int k = beg; k < end; ++k) wls_c[as_size(ap_flat[as_size(k)])] -= delta;
                }
                polished_best = std::min(polished_best, slack_objective());
                if (!changed) break;
            }
            best = polished_best;
        }

        store_best_dual();
        return static_cast<long double>(task.current_gain) + eval_dense_exact(active_best_y, true);
    }

    [[nodiscard]] long double lagrangian_upper_bound(const Task& task, WorkerContext& ctx, int* greedy_gain_out = nullptr) const {
        return lagrangian_upper_bound(task, ctx, 64, greedy_gain_out); // optimal value ranges found around 62-66 and 112-116. 64 was the best single value in testing.
    }

    void capture_root_cost_lb(long double gain_ub) {
        if (root_lb_captured_) return;
        raw_root_cost_lb_ = cost_from_gain(cur_n_, as_int(std::floor(gain_ub + 1e-9L)));
        root_cost_lb_ = std::max(raw_root_cost_lb_, root_cost_lb_floor_);
        root_lb_captured_ = true;
    }

    void maybe_submit_greedy(const Task& task, WorkerContext& ctx, int greedy_gain = -1, bool greedy_lines_ready = false) {
        if (!greedy_lines_ready) greedy_gain = greedy_completion(task, ctx, true);
        if (task.current_gain + greedy_gain <= current_best_gain()) return;

        ctx.candidate_choice = task.current_choice;
        ctx.candidate_choice.insert(ctx.candidate_choice.end(), ctx.greedy_lines.begin(), ctx.greedy_lines.end());
        submit_candidate(task.current_gain + greedy_gain, ctx.candidate_choice);
    }

    [[nodiscard]] long double exclude_branch_upper_bound(const Task& task, const WorkerContext& ctx, int branch_line, long double node_gain_ub) const {
        const auto& best_y = ctx.lagrangian.best_y;
        const auto& line_slack = ctx.lagrangian.line_slack;
        long double capped = node_gain_ub - std::max(0.0L, line_slack[as_size(branch_line)]);

        for (int p : lines_[as_size(branch_line)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            const long double residual = 1.0L - best_y[as_size(p)];
            if (residual <= 1e-18L) continue;

            long double requiring_sum = 0.0L;
            long double requiring_best = 0.0L;
            for (int id : ordered_point_lines_[as_size(p)]) {
                if (id == branch_line || !task.available[as_size(id)] || task.line_cover[as_size(id)] < 3) continue;
                const long double rc = line_slack[as_size(id)];
                if (rc <= 1e-18L || rc > residual + 1e-18L) continue;
                requiring_sum += rc;
                requiring_best = std::max(requiring_best, rc);
            }
            capped -= (requiring_sum - requiring_best);
        }
        return std::max<long double>(static_cast<long double>(task.current_gain), capped);
    }

    void clear_seen_state_table() {
        for (SeenStateShard& shard : seen_state_shards_) shard.best_gain.clear();
    }

    void collect_blocked_productive_lines(const Task& task, std::vector<int>& out) const {
        out = task.blocked_productive_ids;  // already sorted
    }

    [[nodiscard]] StateKey make_state_key(const Task& task, WorkerContext& ctx) const {
        collect_blocked_productive_lines(task, ctx.state_blocked_ids);
        return {task.active, ctx.state_blocked_ids};
    }

    [[nodiscard]] bool seen_search_state(const Task& task, WorkerContext& ctx) {
        StateKey key = make_state_key(task, ctx);
        const std::size_t hash = StateKeyHash{}(key);
        SeenStateShard& shard = seen_state_shards_[hash % kSeenStateShardCount];
        std::scoped_lock lock(shard.mutex);
        if (auto it = shard.best_gain.find(key); it != shard.best_gain.end()) {
            if (it->second >= task.current_gain) return true;
            it->second = task.current_gain;
            return false;
        }
        shard.best_gain.emplace(std::move(key), task.current_gain);
        return false;
    }

    [[nodiscard]] bool remember_frontier_state(const Task& task, WorkerContext& ctx, StateBestMap& best_gain) const {
        StateKey key = make_state_key(task, ctx);
        if (auto it = best_gain.find(key); it != best_gain.end()) {
            if (it->second >= task.current_gain) return true;
            it->second = task.current_gain;
            return false;
        }
        best_gain.emplace(std::move(key), task.current_gain);
        return false;
    }

    [[nodiscard]] bool frontier_state_is_stale(const Task& task, WorkerContext& ctx, const StateBestMap& best_gain) const {
        StateKey key = make_state_key(task, ctx);
        auto it = best_gain.find(key);
        return it != best_gain.end() && it->second > task.current_gain;
    }

    void dfs(Task& task, WorkerContext& ctx, int depth = 0) {
        if (root_interval_closed()) return;
        int observed = depth_max_atomic_ptr->load(std::memory_order_relaxed);
        while (depth > observed &&
               !depth_max_atomic_ptr->compare_exchange_weak(observed, depth, std::memory_order_relaxed, std::memory_order_relaxed)) {}
        const UndoMark node_undo = mark_undo(ctx);
        const int start_gain = task.current_gain;
        const std::size_t start_choice_size = task.current_choice.size();
        
        auto rollback = [&]() {
            task.current_choice.resize(start_choice_size);
            task.current_gain = start_gain;
            undo_to(task, ctx, node_undo);
        };

        apply_forced_line_reduction(task, ctx, true);
        node_count_.fetch_add(1, std::memory_order_relaxed);
        if (seen_search_state(task, ctx)) {
            seen_prune_atomic_ptr->fetch_add(1, std::memory_order_relaxed);
            rollback();
            return;
        }
        
        int greedy_gain = 0;
        const int lag_max_iters = std::min(64, 8 + task.active_count / 8);
        const long double ub = lagrangian_upper_bound(task, ctx, lag_max_iters, &greedy_gain);
        capture_root_cost_lb(ub);
        
        const bool lag_killed = gain_upper_bound_cannot_beat_best_cost(ub);
        if (root_interval_closed() || lag_killed) {
            if (lag_killed) lag_prune_atomic_ptr->fetch_add(1, std::memory_order_relaxed);
            rollback();
            return;
        }

        maybe_submit_greedy(task, ctx, greedy_gain, true);
        if (root_interval_closed()) { rollback(); return; }

        const std::size_t dual_mark = save_current_dual(ctx);
        const BranchChoice branch = choose_branch(task, ctx, ub);
        const int branch_line = branch.line;
        
        if (branch_line < 0) {
            restore_saved_dual(ctx, dual_mark);
            submit_candidate(task.current_gain, task.current_choice);
            rollback(); return;
        }

        const long double exclude_ub = branch.exclude_ub;

        {
            const UndoMark undo = mark_undo(ctx);
            const int delta = positive_gain(task, branch_line);
            if (delta > 0) {
                select_line(task, ctx, branch_line);
                task.current_gain += delta;
                task.current_choice.push_back(branch_line);
                dfs(task, ctx, depth + 1);
                restore_saved_dual(ctx, dual_mark);
                task.current_choice.pop_back();
                task.current_gain -= delta;
            }
            undo_to(task, ctx, undo);
        }

        if (gain_upper_bound_cannot_beat_best_cost(exclude_ub)) { rollback(); return; }

        {
            const UndoMark undo = mark_undo(ctx);
            disable_line(task, ctx, branch_line);
            dfs(task, ctx, depth + 1);
            restore_saved_dual(ctx, dual_mark);
            undo_to(task, ctx, undo);
        }

        rollback();
    }

    void run_task(Task task, WorkerContext& ctx) {
        ensure_task_search_buffers(task);
        prepare_worker_context(ctx);
        if (root_interval_closed()) return;
        apply_forced_line_reduction(task, ctx, false);
        maybe_submit_greedy(task, ctx);
        if (root_interval_closed()) return;
        dfs(task, ctx);
    }

    [[nodiscard]] std::vector<Task> build_frontier(Task root, WorkerContext& ctx) {
        std::vector<Task> frontier;
        frontier.reserve(std::max<std::size_t>(4, as_size(worker_count_) * 8U));
        frontier.push_back(std::move(root));

        std::vector<Task> children;
        const unsigned frontier_multiplier = current_best_cost() >= 93 ? 8U : 4U;
        const std::size_t target = std::max<std::size_t>(1, as_size(worker_count_) * as_size(frontier_multiplier));
        StateBestMap frontier_best;
        frontier_best.reserve(target * 4);
        static_cast<void>(remember_frontier_state(frontier.front(), ctx, frontier_best));

        while (frontier.size() < target) {
            if (root_interval_closed()) { frontier.clear(); break; }
            
            auto it = std::ranges::max_element(frontier, [](const Task& a, const Task& b) {
                if (a.active_count != b.active_count) return a.active_count < b.active_count;
                return a.current_gain > b.current_gain;
            });
            if (it == frontier.end()) break;

            Task task = std::move(*it);
            frontier.erase(it);
            apply_forced_line_reduction(task, ctx, false);
            if (frontier_state_is_stale(task, ctx, frontier_best)) continue;

            ctx.lagrangian.dual_ready = false;
            int greedy_gain = 0;
            const int lag_max_iters = std::min(64, 8 + task.active_count / 8);
            const long double ub = lagrangian_upper_bound(task, ctx, lag_max_iters, &greedy_gain);
            capture_root_cost_lb(ub);
            if (root_interval_closed()) { frontier.clear(); break; }
            if (gain_upper_bound_cannot_beat_best_cost(ub)) {
                if (frontier.empty()) break;
                continue;
            }

            maybe_submit_greedy(task, ctx, greedy_gain, true);
            if (root_interval_closed()) { frontier.clear(); break; }

            const BranchChoice branch = choose_branch(task, ctx, ub);
            const int branch_line = branch.line;
            if (branch_line < 0) {
                submit_candidate(task.current_gain, task.current_choice);
                if (frontier.empty()) break;
                continue;
            }

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

            if (!gain_upper_bound_cannot_beat_best_cost(branch.exclude_ub)) {
                Task exclude = task;
                ensure_task_search_buffers(exclude);
                disable_line_without_undo(exclude, branch_line);
                children.push_back(std::move(exclude));
            }

            for (Task& child : children) {
                apply_forced_line_reduction(child, ctx, false);
                if (!remember_frontier_state(child, ctx, frontier_best)) frontier.push_back(std::move(child));
            }
        }

        std::ranges::sort(frontier, [](const Task& a, const Task& b) {
            if (a.active_count != b.active_count) return a.active_count > b.active_count;
            if (a.productive_ids.size() != b.productive_ids.size()) return a.productive_ids.size() > b.productive_ids.size();
            return a.current_gain < b.current_gain;
        });

        return frontier;
    }

    void solve_frontier_parallel(const std::vector<Task>& frontier) {
        const unsigned threads = std::min<unsigned>(worker_count_, static_cast<unsigned>(frontier.size()));
        if (threads <= 1) { WorkerContext ctx; run_task(frontier.front(), ctx); return; }

        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        workers.reserve(as_size(threads - 1));

        auto worker = [&]() {
            WorkerContext ctx;
            while (true) {
                const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= frontier.size()) break;
                run_task(frontier[idx], ctx);
            }
        };

        for (unsigned i = 1; i < threads; ++i) workers.emplace_back(worker);
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
    bool                                    trace_this_n_ = false;
    bool                                    force_exact_forced_lines_ = false;
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
}; // end ExactGainSolver

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
struct CoverInequalityResult {
    int                        lb = 0;
    std::vector<unsigned char> in_support;
};

[[nodiscard]] CoverInequalityResult cover_inequality_lb(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
{
    CoverInequalityResult result;
    if (n == 0) return result;

    result.in_support.assign(as_size(n), 1);
    if (active_ids.empty()) {
        result.lb = (n + 1) / 2;
        return result;
    }

    const int m = as_int(active_ids.size());
    // local_inc[p] = indices into active_ids for lines containing point p.
    std::vector<std::vector<int>> local_inc(as_size(n));
    std::vector<int> line_size(as_size(m), 0);
    
    for (int i = 0; i < m; ++i) {
        for (int p : lines[as_size(active_ids[i])].points) {
            if (p >= n) break;
            local_inc[as_size(p)].push_back(i);
            ++line_size[i];
        }
    }

    // S starts as all n points.  violation_count[p] = number of active lines
    // through p that currently have |L ∩ S| >= 3.
    int S_size = n;
    std::vector<int> violation_count(as_size(n), 0);
    for (int i = 0; i < m; ++i) {
        if (line_size[i] >= 3) {
            for (int p : lines[as_size(active_ids[i])].points) {
                if (p >= n) break;
                ++violation_count[as_size(p)];
            }
        }
    }

    // Greedy: remove the point with the most violated-line incidences.
    while (true) {
        int worst_p = -1, worst_v = 0;
        for (int p = 0; p < n; ++p) {
            if (result.in_support[as_size(p)] && violation_count[as_size(p)] > worst_v) {
                worst_v = violation_count[as_size(p)];
                worst_p = p;
            }
        }
        if (worst_v == 0) break;

        result.in_support[as_size(worst_p)] = 0;
        --S_size;

        for (int i : local_inc[as_size(worst_p)]) {
            const bool was_violated = (line_size[i] >= 3);
            if (--line_size[i] < 3 && was_violated) {
                // Line i just fell below the threshold; un-charge its S members.
                for (int q : lines[as_size(active_ids[i])].points) {
                    if (q >= n) break;
                    if (result.in_support[as_size(q)]) --violation_count[as_size(q)];
                }
            }
        }
    }

    result.lb = (S_size + 1) / 2;
    return result;
}

[[nodiscard]] ExactSolveResult solve_exact(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines,
    const std::vector<std::vector<int>>& incidence, int root_lb_floor,
    const std::vector<int>& greedy_seed_line_ids, std::vector<int>& warm_start_solution)
{
    // Strengthen the root lower-bound floor with the cover inequality before
    // handing off to the gain solver. If this closes the gap at the root,
    // root_interval_closed() fires after the first Lagrangian evaluation and
    // the DFS is eliminated entirely with zero branching.
    ExactSolveResult out;
    const CoverInequalityResult cover = cover_inequality_lb(n, active_ids, lines);
    const int effective_lb_floor = std::max(root_lb_floor, cover.lb);
    ExactGainSolver solver(lines, incidence);
    GainSolveResult result = solver.solve_for_n(n, active_ids, warm_start_solution, greedy_seed_line_ids, effective_lb_floor);
    
    out.answer = (n - result.best_gain + 1) / 2;
    out.lb = result.root_cost_lb;
    out.lb_raw = result.raw_root_cost_lb;
    out.lb_cov = cover.lb;
    out.lb_lp = 0;
    out.prod = result.root_productive_lines;
    out.pinc = result.root_productive_incidence;
    out.nodes = result.nodes;
    out.chosen_line_ids = std::move(result.chosen_line_ids);
    out.timed_out = result.timed_out;
    warm_start_solution = out.chosen_line_ids;
    return out;
}

} // namespace solver

// ---------------------------------------------------------------------------
// User Interface
// ---------------------------------------------------------------------------
namespace ui {
    using util::as_size;
    using util::as_int;

// ---------------------------------------------------------------------------
// LiveDisplay — background thread that prints running progress to stderr.
// Only activates after a single N has been solving for more than 0.2 second,
// so it never clutters output for fast primes. All writes go to stderr;
// stdout carries only the final per-N result lines.
//
// The main thread calls begin_n() just before each solve and clear() just
// after timing is captured. begin_n() remains two relaxed atomic stores on the
// hot path, while clear() runs outside the timed region and also resets the
// active N so the display thread cannot redraw a stale status line.
// ---------------------------------------------------------------------------
class LiveDisplay {
public:
    using Clock = std::chrono::steady_clock;

    explicit LiveDisplay(Clock::time_point program_start)
        : program_start_(program_start), thread_([this] { run(); }) {}

    ~LiveDisplay() {
        active_.store(false, std::memory_order_relaxed);
        if (thread_.joinable()) thread_.join();
        clear_line();
    }

    void begin_n(int n) noexcept {
        n_start_ns_.store(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - program_start_).count(), std::memory_order_relaxed);
        current_n_.store(n, std::memory_order_relaxed);
    }

    void clear() const noexcept {
        current_n_.store(0, std::memory_order_relaxed);
        n_start_ns_.store(0, std::memory_order_relaxed);
        clear_line();
    }

private:
    static constexpr double kActivationDelaySeconds = 0.2;

    void clear_line() const {
        std::print(stderr, "\r                                                                              \r");
        std::fflush(stderr);
    }

    void run() {
        while (active_.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(100ms);
            if (!active_.load(std::memory_order_relaxed)) break;

            const int n = current_n_.load(std::memory_order_relaxed);
            if (n == 0) continue;

            const std::int64_t start_ns = n_start_ns_.load(std::memory_order_relaxed);
            if (start_ns == 0) continue;

            const std::int64_t now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - program_start_).count();
            const double n_elapsed = static_cast<double>(now_ns - start_ns) * 1e-9;
            if (n_elapsed < kActivationDelaySeconds) continue;

            std::print(stderr, "\r[running] N={}  N_time={:.1f}s  total={:.1f}s  ", n, n_elapsed, static_cast<double>(now_ns) * 1e-9);
            std::fflush(stderr);
        }
    }

    Clock::time_point                 program_start_;
    std::atomic<bool>                 active_{true};
    mutable std::atomic<int>          current_n_{0};
    mutable std::atomic<std::int64_t> n_start_ns_{0};
    std::thread                       thread_;
};

} // namespace ui
} // namespace

int main(int argc, char** argv) {
    using namespace config;
    using namespace util;
    using namespace debug;
    using namespace core;
    using namespace solver;
    using namespace ui;

    const auto total_start = std::chrono::steady_clock::now();
    int requested_n = kExecutionLimit;
    if (argc >= 2) {
        if (const int parsed = std::atoi(argv[1]); parsed > 0) requested_n = parsed;
    }
    requested_n = std::min({requested_n, kExecutionLimit, kBitCapacity});

    const std::vector<int>       primes = generate_primes(requested_n);
    const std::vector<HeavyLine> lines  = enumerate_heavy_lines(primes);

    std::vector<std::vector<int>> incidence(as_size(requested_n));
    std::vector<std::vector<int>> activation(as_size(requested_n + 1));
    for (int id = 0; id < as_int(lines.size()); ++id) {
        for (int p : lines[as_size(id)].points) {
            if (p < requested_n) incidence[as_size(p)].push_back(id);
        }
        if (lines[as_size(id)].activate_at <= requested_n) {
            activation[as_size(lines[as_size(id)].activate_at)].push_back(id);
        }
    }

    std::vector<int> active_ids;
    active_ids.reserve(lines.size());

    const auto program_start = std::chrono::steady_clock::now();
    LiveDisplay display(program_start);

    int  prev = (kStartN > 1) ? kStartN - 1 : 0;
    long long total_nodes_processed = 0;
    std::atomic<int> forced_atomic{0};
    std::atomic<long long> lag_iters_atomic{0};
    std::atomic<int> depth_max_atomic{0};
    std::atomic<int> seen_prune_atomic{0};
    std::atomic<int> lag_prune_atomic{0};
    std::atomic<int> strong_branch_atomic{0};
    std::vector<int> warm_start_solution;
    WitnessCover witness_cover;
    bool witness_ready = false;

    forced_atomic_ptr = &forced_atomic;
    lag_iters_atomic_ptr = &lag_iters_atomic;
    depth_max_atomic_ptr = &depth_max_atomic;
    seen_prune_atomic_ptr = &seen_prune_atomic;
    lag_prune_atomic_ptr = &lag_prune_atomic;
    strong_branch_atomic_ptr = &strong_branch_atomic;
    
    for (int n_init = 1; n_init < kStartN; ++n_init) {
        for (int id : activation[as_size(n_init)]) active_ids.push_back(id);
    }

    for (int n = kStartN; n <= requested_n; ++n) {
        for (int id : activation[as_size(n)]) active_ids.push_back(id);
        forced_atomic.store(0, std::memory_order_relaxed);
        lag_iters_atomic.store(0, std::memory_order_relaxed);
        depth_max_atomic.store(0, std::memory_order_relaxed);
        seen_prune_atomic.store(0, std::memory_order_relaxed);
        lag_prune_atomic.store(0, std::memory_order_relaxed);
        strong_branch_atomic.store(0, std::memory_order_relaxed);

        display.begin_n(n);
        const auto t0 = std::chrono::steady_clock::now();

        const GreedyCoverSolution greedy = greedy_upper_bound_solution(n, active_ids, lines);
        const int ub_raw = greedy.cost;
        int incumbent = std::min({n, prev + 1, greedy.cost});
        const int lb0_floor = warm_start_solution.empty() ? 0 : prev;
        const std::size_t active_count = active_ids.size();
        const bool witness_hit = witness_ready && witness_covers_point(witness_cover, n - 1, primes[as_size(n - 1)]);
        if (witness_hit) incumbent = std::min(incumbent, prev);
        const int ub0_val = incumbent;

        ExactSolveResult exact_result;
        if (witness_hit) {
            exact_result.answer = prev;
            exact_result.lb = prev;
            exact_result.lb_raw = prev;
        } else {
            exact_result = solve_exact(n, active_ids, lines, incidence, lb0_floor, greedy.line_ids, warm_start_solution);
            if (!exact_result.timed_out) {
                witness_cover = build_witness_cover(n, primes, lines, warm_start_solution);
                witness_ready = true;
            }
        }
        
        const int gap = ub0_val - exact_result.lb;
        const int gap_raw = ub_raw - exact_result.lb_raw;
        const char mode = witness_hit ? 'W' : (exact_result.nodes == 0 ? 'R' : 'D');
        const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        display.clear();
        total_nodes_processed += exact_result.nodes;

        const int forced = forced_atomic.load(std::memory_order_relaxed);
        const long long lag_iters = mode == 'D' ? lag_iters_atomic.load(std::memory_order_relaxed) : 0;
        const int depth = mode == 'D' ? depth_max_atomic.load(std::memory_order_relaxed) : 0;
        const int seen_prune = mode == 'D' ? seen_prune_atomic.load(std::memory_order_relaxed) : 0;
        const int lag_prune = mode == 'D' ? lag_prune_atomic.load(std::memory_order_relaxed) : 0;
        const int strong_branch = mode == 'D' ? strong_branch_atomic.load(std::memory_order_relaxed) : 0;

        if (exact_result.timed_out) {
            std::println("[stopped] N={} time limit of {}s exceeded (cancelled mid-run after {})", 
                         n, kPerNTimeLimitSeconds, format_seconds(elapsed));
            break;
        }

        if (debug_trace_enabled_for_n(n)) {
            const ChoiceReplay replay = replay_choice_trace(n, lines, exact_result.chosen_line_ids);
            debug_write_text(".chosen.txt", format_choice_trace(n, exact_result.answer, replay, exact_result.chosen_line_ids), "w");
        }

        std::println("N={} prime={} active={} mode={} lb_floor={} lb_cov={} lb_lp={}"
                     " ub0={} lb={} gap={} ub_raw={} lb_raw={} gap_raw={} forced={}"
                     " prod={} pinc={} lag_iters={} depth={} seen_prune={} lag_prune={} strong_branch={}"
                     " lines={} time={} nodes={}",
            n,
            primes[as_size(n - 1)],
            active_count,
            mode,
            lb0_floor,
            exact_result.lb_cov,
            exact_result.lb_lp,
            ub0_val,
            exact_result.lb,
            gap,
            ub_raw,
            exact_result.lb_raw,
            gap_raw,
            forced,
            exact_result.prod,
            exact_result.pinc,
            lag_iters,
            depth,
            seen_prune,
            lag_prune,
            strong_branch,
            exact_result.answer,
            format_seconds(elapsed),
            exact_result.nodes);
        
        prev = exact_result.answer;

        if (kPerNTimeLimitSeconds > 0.0 && elapsed > kPerNTimeLimitSeconds) {
            std::println("[stopped] per-N time limit of {}s exceeded", kPerNTimeLimitSeconds);
            break;
        }
    }
    const double total_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start).count();
    std::println("total time={} total nodes={}", format_seconds(total_elapsed), total_nodes_processed);
    return 0;
}
