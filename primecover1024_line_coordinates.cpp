/*
 * primecover1024.cpp
 *
 * Exact minimum line cover solver for prime points (i, p_i), 1 <= i <= N.
 *
 * Features:
 * - Enumerates heavy lines (lines containing >= 3 points).
 * - Solves residual exact-gain search using 1024-bit bitmasks.
 * - Applies greedy seeding, a Lagrangian upper bound, and a shared atomic incumbent.
 * - Uses a root cover-inequality lower bound before branch-and-bound.
 * - Closes uncovered residues with 1- or 2-point lines at cost ceil(residual / 2).
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
#include <limits>
#include <mutex>
#include <numeric>
#include <print>
#include <ranges>
#include <span>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#ifndef PC_SKIP_FORCE_REBUILD_PRODUCTIVE_METADATA
#define PC_SKIP_FORCE_REBUILD_PRODUCTIVE_METADATA 1
#endif

#ifndef PC_VALIDATE_FORCE_METADATA
#define PC_VALIDATE_FORCE_METADATA 0
#endif

namespace {

// ===========================================================================
// Configuration
// ===========================================================================
namespace config {
    constexpr int    kBitCapacity          = 1024; // DO NOT CHANGE - the core solver relies on this for fixed-size bitmask operations.
    constexpr int    kBitWords             = kBitCapacity / 64; // DO NOT CHANGE - derived from kBitCapacity.
    constexpr int    kStartN               = 0; // Move freely - Starting N for the sweep; 0 to start from 1
    constexpr int    kExecutionLimit       = 1024; // Move freely - the solver stops when it reaches this N. (1024 limit)
    
    // Wall-clock timeout limit per N in seconds. Set to 0 for no limit.
    constexpr double kPerNTimeLimitSeconds = 0; // Move freely
} // namespace config

// ===========================================================================
// Utility Helpers
// ===========================================================================
namespace util {
    template <typename T>
    [[nodiscard]] constexpr std::size_t as_size(T v) noexcept {
        return static_cast<std::size_t>(v);
    }

    template <typename T>
    [[nodiscard]] constexpr int as_int(T v) noexcept {
        return static_cast<int>(v);
    }

    [[nodiscard]] constexpr int words_for_n(int n) noexcept { return (n + 63) / 64; }
    [[nodiscard]] constexpr int ceil_half(int v)   noexcept { return (v + 1) / 2; }

    [[nodiscard]] inline std::string format_seconds(double s) {
        return std::format("{:.{}f}s", s, s < 10.0 ? 6 : 3);
    }

    [[nodiscard]] constexpr std::uint64_t mix_u64(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x  = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x  = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    [[nodiscard]] constexpr std::int64_t abs64(std::int64_t x) noexcept { 
        return x < 0 ? -x : x; 
    }
} // namespace util

// ===========================================================================
// Geometry
// ===========================================================================
namespace geom {
    using util::as_size;
    using util::mix_u64;
    using util::abs64;

    struct LineKey {
        std::int64_t a{0}, b{0}, c{0};
        bool operator==(const LineKey&) const = default;
    };

    struct SlopeKey {
        int dy{0};
        int dx{1};
        bool operator==(const SlopeKey&) const = default;
    };

    struct SlopeKeyHash {
        [[nodiscard]] std::size_t operator()(const SlopeKey& k) const noexcept {
            const std::uint64_t packed =
                (static_cast<std::uint64_t>(static_cast<std::uint32_t>(k.dy)) << 32) |
                static_cast<std::uint32_t>(k.dx);
            return as_size(mix_u64(packed));
        }
    };

    // Normalizes line equation ax + by + c = 0
    [[nodiscard]] LineKey canonical_line(int x1, int y1, int x2, int y2) noexcept {
        std::int64_t a = static_cast<std::int64_t>(y1) - static_cast<std::int64_t>(y2);
        std::int64_t b = static_cast<std::int64_t>(x2) - static_cast<std::int64_t>(x1);
        std::int64_t c = -(a * static_cast<std::int64_t>(x1) + b * static_cast<std::int64_t>(y1));
        
        if (const std::int64_t g = std::gcd(std::gcd(abs64(a), abs64(b)), abs64(c)); g) { 
            a /= g; b /= g; c /= g; 
        }
        
        if (a < 0 || (a == 0 && b < 0)) {
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

// ===========================================================================
// Core Data Structures & Generation
// ===========================================================================
namespace core {
    using namespace geom;
    using util::as_size;
    using util::as_int;
    using util::ceil_half;
    using util::words_for_n;

    struct BitMask1024 {
        std::array<std::uint64_t, config::kBitWords> words{};

        void set(int bit) noexcept   { words[as_size(bit >> 6)] |=  (1ULL << (bit & 63)); }
        void reset(int bit) noexcept { words[as_size(bit >> 6)] &= ~(1ULL << (bit & 63)); }
        
        [[nodiscard]] bool test(int bit) const noexcept {
            return ((words[as_size(bit >> 6)] >> (bit & 63)) & 1ULL) != 0;
        }

        [[nodiscard]] int intersection_count(const BitMask1024& other, int used_words) const noexcept {
            int total = 0;
            for (int i = 0; i < used_words; ++i) {
                total += std::popcount(words[as_size(i)] & other.words[as_size(i)]);
            }
            return total;
        }

        [[nodiscard]] BitMask1024 and_not(const BitMask1024& other) const noexcept {
            BitMask1024 out;
            for (int i = 0; i < config::kBitWords; ++i) {
                out.words[as_size(i)] = words[as_size(i)] & ~other.words[as_size(i)];
            }
            return out;
        }

        [[nodiscard]] static BitMask1024 prefix(int n) noexcept {
            BitMask1024 out;
            const int full_words = n / 64;
            for (int i = 0; i < full_words; ++i) {
                out.words[as_size(i)] = std::numeric_limits<std::uint64_t>::max();
            }
            if (const int tail = n % 64; tail != 0) {
                out.words[as_size(full_words)] = (1ULL << tail) - 1ULL;
            }
            return out;
        }

        template <class Fn>
        void for_each_set_bit(int used_words, Fn&& fn) const {
            for (int wi = 0; wi < used_words; ++wi) {
                std::uint64_t word = words[as_size(wi)];
                while (word != 0) {
                    fn((wi << 6) + std::countr_zero(word));
                    word &= word - 1;
                }
            }
        }
    };

    struct HeavyLine {
        BitMask1024      mask{};
        std::vector<int> points;
        int              activate_at{0};
    };

    struct GreedyCoverSolution {
        int              cost{0};
        std::vector<int> line_ids;
    };

    struct GainSolveResult {
        int              best_gain{0};
        int              root_cost_lb{0};
        int              raw_root_cost_lb{0};
        int              root_productive_lines{0}; 
        int              root_productive_incidence{0}; 
        int              forced{0};
        // Lines forced by the exclusive dependency rule in the initial root preprocessing pass only, before frontier build or DFS.
        int              forced_root{0};
        long long        lag_iters{0};
        // Total coordinate descent polish sweeps across all Lagrangian evaluations in this solve.
        long long        polish_sweeps{0};
        int              depth_max{0};
        int              lag_prune{0};
        int              strong_branch{0};
        long long        nodes{0};
        // Number of independent subtrees sent to parallel workers; 1 means the single-thread fallback ran.
        int              frontier_size{0};
        // Fractional part of the root Lagrangian gain bound; near 1.0 means lb_raw is close to its next tighter integer.
        double           root_ub_frac{0.0};
        std::vector<int> chosen_line_ids;
        std::vector<double> dual_seed;
        bool             timed_out{false};
    };

    struct ExactSolveResult {
        int              answer{0};
        int              lb{0};
        int              lb_raw{0};
        int              lb_cov{0}; // Root lower bound from cover_inequality_lb()
        int              prod{0}; 
        int              pinc{0}; 
        int              forced{0};
        // Lines forced by the exclusive dependency rule in the initial root preprocessing pass only, before frontier build or DFS.
        int              forced_root{0};
        long long        lag_iters{0};
        // Total coordinate descent polish sweeps across all Lagrangian evaluations in this solve.
        long long        polish_sweeps{0};
        int              depth{0};
        int              lag_prune{0};
        int              strong_branch{0};
        long long        nodes{0};
        // Number of independent subtrees sent to parallel workers; 1 means the single-thread fallback ran.
        int              frontier_size{0};
        // Number of heavy lines from warm_start_solution carried in from the previous N; zero means cold start.
        int              warm_size{0};
        // Heavy lines (>=3 points) selected by greedy_upper_bound_solution; ub_raw minus this is 2-pt/1-pt padding cost.
        int              greedy_heavy{0};
        // Fractional part of the root Lagrangian gain bound; near 1.0 means lb_raw is close to its next tighter integer.
        double           root_ub_frac{0.0};
        std::vector<int> chosen_line_ids;
        bool             timed_out{false};
    };

    struct WitnessCover {
        std::vector<LineKey> lines;
    };

    // Sieve to generate primes up to an estimated maximum bound based on 'count'
    [[nodiscard]] std::vector<int> generate_primes(int count) {
        if (count <= 0) return {};
        
        const double n = static_cast<double>(count);
        const int limit = (count < 10)
            ? 64
            : as_int(n * (std::log(n) + std::log(std::log(n)) + 3.0)) + 256;
            
        std::vector<bool> is_p(as_size(limit + 1), true);
        is_p[0] = is_p[1] = false;
        
        for (int p = 2; p * p <= limit; ++p) {
            if (is_p[as_size(p)]) {
                for (int q = p * p; q <= limit; q += p) {
                    is_p[as_size(q)] = false;
                }
            }
        }
                    
        std::vector<int> primes;
        primes.reserve(as_size(count));
        for (int v = 2; v <= limit && as_int(primes.size()) < count; ++v) {
            if (is_p[as_size(v)]) primes.push_back(v);
        }
            
        return primes;
    }

    // Ensures we don't build duplicate lines by enforcing that 'anchor' is the origin point.
    [[nodiscard]] bool anchor_is_leftmost_for_slope(
        const std::vector<int>& primes, int anchor, const SlopeKey& slope) noexcept
    {
        const int step_x = slope.dx;
        const std::int64_t step_y = slope.dy;
        std::int64_t expected_y = static_cast<std::int64_t>(primes[as_size(anchor)]) - step_y;
        
        for (int x = anchor - step_x; x >= 0; x -= step_x, expected_y -= step_y) {
            if (expected_y == static_cast<std::int64_t>(primes[as_size(x)])) {
                return false;
            }
        }
        return true;
    }

    // Identifies all valid lines containing at least 3 points.
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

    // Generates a greedy cover. Repeatedly selects the line covering the most active points.
    [[nodiscard]] GreedyCoverSolution greedy_upper_bound_solution(
        int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
    {
        const int words = words_for_n(n);
        BitMask1024 active = BitMask1024::prefix(n);
        GreedyCoverSolution result;
        int count = n;
        
        while (true) {
            int best = -1;
            int bc = 2;
            
            for (int id : active_ids) {
                const int c = lines[as_size(id)].mask.intersection_count(active, words);
                if (c > bc || (c == bc && best != -1 && lines[as_size(id)].points.size() > lines[as_size(best)].points.size())) { 
                    bc = c; 
                    best = id; 
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
        int n, int future_end, const std::vector<int>& primes, const std::vector<HeavyLine>& lines, const std::vector<int>& chosen_heavy_line_ids)
    {
        WitnessCover witness;
        witness.lines.reserve(chosen_heavy_line_ids.size() + as_size(ceil_half(n)));

        BitMask1024 covered;
        for (int id : chosen_heavy_line_ids) {
            const auto& heavy = lines[as_size(id)];
            witness.lines.push_back(canonical_line(
                heavy.points[0] + 1, primes[as_size(heavy.points[0])],
                heavy.points[1] + 1, primes[as_size(heavy.points[1])]));
                
            for (int p : heavy.points) {
                if (p >= n) break;
                covered.set(p);
            }
        }

        std::vector<int> residual_points;
        residual_points.reserve(as_size(n));
        for (int p = 0; p < n; ++p) {
            if (!covered.test(p)) residual_points.push_back(p);
        }

        if (future_end > n && residual_points.size() >= 2) {
            std::vector<unsigned char> available(residual_points.size(), 1);
            std::vector<unsigned char> future_covered(as_size(future_end + 1), 0);

            auto update_future_covered = [&](const LineKey& line) {
                for (int t = n + 1; t <= future_end; ++t) {
                    if (!future_covered[as_size(t)] &&
                        line_contains_point(line, t, primes[as_size(t - 1)])) {
                        future_covered[as_size(t)] = 1;
                    }
                }
            };

            for (const LineKey& line : witness.lines) {
                update_future_covered(line);
            }

            while (true) {
                int next_miss = -1;
                for (int t = n + 1; t <= future_end; ++t) {
                    if (!future_covered[as_size(t)]) {
                        next_miss = t;
                        break;
                    }
                }
                if (next_miss < 0) break;

                int best_i = -1;
                int best_j = -1;
                int best_reach = -1;
                int best_hits = -1;
                LineKey best_line{};

                for (std::size_t i = 0; i < residual_points.size(); ++i) {
                    if (!available[i]) continue;
                    for (std::size_t j = i + 1; j < residual_points.size(); ++j) {
                        if (!available[j]) continue;

                        const int p_i = residual_points[i];
                        const int p_j = residual_points[j];
                        const LineKey candidate = canonical_line(
                            p_i + 1, primes[as_size(p_i)],
                            p_j + 1, primes[as_size(p_j)]);

                        if (!line_contains_point(candidate, next_miss, primes[as_size(next_miss - 1)])) {
                            continue;
                        }

                        int reach = 0;
                        for (int t = next_miss; t <= future_end; ++t) {
                            if (future_covered[as_size(t)] ||
                                line_contains_point(candidate, t, primes[as_size(t - 1)])) {
                                ++reach;
                            } else {
                                break;
                            }
                        }

                        int hits = 0;
                        for (int t = next_miss; t <= future_end; ++t) {
                            if (!future_covered[as_size(t)] &&
                                line_contains_point(candidate, t, primes[as_size(t - 1)])) {
                                ++hits;
                            }
                        }

                        const bool better = (reach > best_reach) ||
                                            (reach == best_reach && hits > best_hits) ||
                                            (reach == best_reach && hits == best_hits &&
                                             (best_i < 0 || p_i < residual_points[as_size(best_i)] ||
                                              (p_i == residual_points[as_size(best_i)] &&
                                               p_j < residual_points[as_size(best_j)])));
                        if (better) {
                            best_i = as_int(i);
                            best_j = as_int(j);
                            best_reach = reach;
                            best_hits = hits;
                            best_line = candidate;
                        }
                    }
                }

                if (best_i < 0) break;

                available[as_size(best_i)] = 0;
                available[as_size(best_j)] = 0;
                witness.lines.push_back(best_line);
                update_future_covered(best_line);
            }

            int pending = -1;
            for (std::size_t i = 0; i < residual_points.size(); ++i) {
                if (!available[i]) continue;
                const int p = residual_points[i];
                if (pending < 0) {
                    pending = p;
                    continue;
                }
                witness.lines.push_back(canonical_line(
                    pending + 1, primes[as_size(pending)],
                    p + 1,       primes[as_size(p)]));
                pending = -1;
            }

            if (pending >= 0) {
                witness.lines.push_back({1, 0, -static_cast<std::int64_t>(pending + 1)});
            }
            return witness;
        }

        int pending = -1;
        for (const int p : residual_points) {
            if (pending < 0) {
                pending = p;
                continue;
            }
            witness.lines.push_back(canonical_line(
                pending + 1, primes[as_size(pending)],
                p + 1,       primes[as_size(p)]));
            pending = -1;
        }

        if (pending >= 0) {
            witness.lines.push_back({1, 0, -static_cast<std::int64_t>(pending + 1)});
        }

        return witness;
    }

    [[nodiscard]] bool witness_covers_point(const WitnessCover& witness, int x, int y) noexcept {
        for (const LineKey& line : witness.lines) {
            if (line_contains_point(line, x, y)) return true;
        }
        return false;
    }

} // namespace core

// ===========================================================================
// Solver Engine
// ===========================================================================
namespace solver {
    using namespace core;
    using util::as_size;
    using util::as_int;
    using util::words_for_n;
    using util::mix_u64;

class ExactGainSolver {
public:
    ExactGainSolver(const std::vector<HeavyLine>& lines, const std::vector<std::vector<int>>& incidence);

    [[nodiscard]] GainSolveResult solve_for_n(
        int n, const std::vector<int>& active_ids, const std::vector<int>& warm_start_line_ids,
        const std::vector<double>& warm_dual_seed,
        const std::vector<int>& greedy_seed_line_ids, int root_cost_lb_floor);

private:
    // -----------------------------------------------------------------------
    // Types & Structs
    // -----------------------------------------------------------------------
    using DualValue         = double;
    using SubgradientValue  = std::int32_t;
    using CompactPointIndex = std::uint16_t;

    static_assert(config::kBitCapacity <= std::numeric_limits<CompactPointIndex>::max());

    struct CompactLine3 { CompactPointIndex a, b, c; };
    struct CompactLine4 { CompactPointIndex a, b, c, d; };
    struct CompactLine5 { CompactPointIndex a, b, c, d, e; };

    struct SearchStats {
        int       forced{0};
        long long lag_iters{0};
        // Total coordinate descent polish sweeps across all Lagrangian evaluations in this solve.
        long long polish_sweeps{0};
        int       depth_max{0};
        int       lag_prune{0};
        int       strong_branch{0};
        long long nodes{0};

        void merge_from(const SearchStats& other) noexcept {
            forced        += other.forced;
            lag_iters     += other.lag_iters;
            polish_sweeps += other.polish_sweeps;
            depth_max      = std::max(depth_max, other.depth_max);
            lag_prune     += other.lag_prune;
            strong_branch += other.strong_branch;
            nodes         += other.nodes;
        }
    };

    struct Task {
        BitMask1024                   active{};
        int                           active_count{0};
        int                           current_gain{0};
        
        std::vector<std::uint16_t>    line_cover;
        std::vector<int>              productive_degree;
        std::vector<int>              sole_productive_line;
        std::vector<unsigned char>    available;
        std::vector<int>              current_choice;
        
        // Compact set of productive line IDs: available[id]==1 && line_cover[id]>=3.
        std::vector<int>              productive_ids;
        std::vector<int>              productive_pos;
        
        // Compact set of unavailable lines with positive residual coverage.
        // Used for duplicate-state pruning.
        std::vector<int>              blocked_productive_ids;
        std::vector<int>              blocked_productive_pos;
    };

    // Snapshots context state for rolling back DFS branches.
    struct UndoMark {
        std::size_t disabled_lines{0};
        std::size_t decremented_lines{0};
        std::size_t newly_covered_points{0};
        std::size_t covered_point_metadata{0};
        std::size_t removed_from_productive{0};
        std::size_t removed_from_blocked_productive{0};
        std::size_t point_metadata_points{0};
        std::size_t undo_token_depth_before{0};
    };

    struct StateKey {
        BitMask1024      active;
        std::vector<int> blocked_productive_ids;

        [[nodiscard]] bool operator==(const StateKey& other) const noexcept {
            return active.words == other.active.words && 
                   blocked_productive_ids == other.blocked_productive_ids;
        }
    };

    struct StateKeyHash {
        [[nodiscard]] std::size_t operator()(const StateKey& key) const noexcept {
            std::uint64_t h = 0x9e3779b97f4a7c15ULL;
            for (const std::uint64_t word : key.active.words) {
                h ^= mix_u64(word + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            }
            for (const int id : key.blocked_productive_ids) {
                h ^= mix_u64(static_cast<std::uint64_t>(id) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2));
            }
            return as_size(h);
        }
    };

    using StateBestMap = std::unordered_map<StateKey, int, StateKeyHash>;

    struct LagrangianScratch {
        std::vector<int>               active_points;
        std::vector<DualValue>         active_y;
        std::vector<DualValue>         active_best_y;
        std::vector<DualValue>         best_y;
        std::vector<SubgradientValue>  active_subgrad;
        
        std::vector<int>               compact_line_ids;
        std::vector<std::size_t>       compact_line_offsets;
        std::vector<CompactPointIndex> compact_line_points;
        std::vector<CompactPointIndex> compact_point_aps;
        std::vector<CompactLine3>      compact_eval3;
        std::vector<CompactLine4>      compact_eval4;
        std::vector<CompactLine5>      compact_eval5;
        std::vector<int>               compact_eval_other_lines;
        bool                           dual_ready{false};

        // Inverse mapping: active-point index (ap) -> compact productive line indices (li).
        // Speeds up coordinate polish by avoiding global incidence scans.
        std::vector<int>               compact_pt_prod_offsets;  
        std::vector<int>               compact_pt_prod_flat;     
        std::vector<int>               compact_ap_fill_pos;      
        std::vector<int>               point_to_ap_scratch;      

        // Compact-indexed slack array for faster vectorization in the polish hot path.
        std::vector<DualValue>         work_line_slack_compact;

        // Precomputed `line_cover[id] - 2` for compact productive lines.
        // Replaces random-access reads into `line_cover` with sequential reads, 
        // reducing L1 cache pressure during tight bound loops.
        std::vector<int>               compact_rem2;             

        void ensure(std::size_t line_count, int current_n) {
            const std::size_t n = as_size(current_n);
            if (active_points.capacity() < n) active_points.reserve(n);
            if (active_y.size() != n) active_y.resize(n);
            if (active_best_y.size() != n) active_best_y.resize(n);
            if (best_y.size() != n) best_y.resize(n);
            if (active_subgrad.size() != n) active_subgrad.resize(n);
            
            if (compact_line_ids.capacity() < line_count) compact_line_ids.reserve(line_count);
            if (compact_line_offsets.capacity() < line_count + 1) compact_line_offsets.reserve(line_count + 1);
            if (compact_eval3.capacity() < line_count) compact_eval3.reserve(line_count);
            if (compact_eval4.capacity() < line_count) compact_eval4.reserve(line_count);
            if (compact_eval5.capacity() < line_count) compact_eval5.reserve(line_count);
            if (compact_eval_other_lines.capacity() < line_count) compact_eval_other_lines.reserve(line_count);
            
            if (as_int(point_to_ap_scratch.size()) < current_n) {
                point_to_ap_scratch.assign(as_size(current_n), -1);
            }
        }
    };

    struct WorkerContext {
        LagrangianScratch                         lagrangian;
        std::array<int, config::kBitCapacity + 1> greedy_bucket_head{};
        std::vector<int>                          greedy_bucket_next;
        std::vector<int>                          greedy_lines;
        std::vector<int>                          candidate_choice;
        std::vector<int>                          disabled_lines;
        std::vector<int>                          decremented_lines;
        std::vector<int>                          newly_covered_points;
        std::vector<int>                          removed_from_productive;
        std::vector<int>                          removed_from_blocked_productive;
        std::vector<int>                          state_blocked_ids;
        std::vector<int>                          saved_dual_points;
        std::vector<DualValue>                    saved_dual_values;
        std::vector<int>                          forced_touched_lines;
        std::vector<int>                          covered_point_degree;
        std::vector<int>                          covered_point_sole;
        std::vector<int>                          include_ub_touched_lines;
        std::vector<DualValue>                    include_ub_line_reduction;
        std::vector<int>                          undo_token_stack;
        int                                       next_undo_token{0};
        std::vector<int>                          point_metadata_points;
        std::vector<int>                          point_metadata_degree;
        std::vector<int>                          point_metadata_sole;
        std::vector<int>                          point_metadata_prev_token;
        std::vector<int>                          point_metadata_saved_token;
        SearchStats                               stats;
    };

    static constexpr int kNoProductiveLine      = -1;
    static constexpr int kManyProductiveLines   = -2;
    static constexpr int kStrongBranchShortlist = 4;
    static constexpr int kStrongBranchMaxIters  = 24;
    static constexpr int kStrongBranchMaxDepth  = 6;

    struct BranchChoice {
        int         line{-1};
        long double exclude_ub{0.0L};
    };

    // -----------------------------------------------------------------------
    // Inline Queries & Utilities
    // -----------------------------------------------------------------------
    [[nodiscard]] int current_best_gain() const noexcept { return best_gain_.load(std::memory_order_relaxed); }
    [[nodiscard]] int current_best_cost() const noexcept { return cost_from_gain(cur_n_, current_best_gain()); }
    [[nodiscard]] bool root_interval_closed() const noexcept { return root_lb_captured_ && current_best_cost() <= root_cost_lb_; }
    [[nodiscard]] int gain_needed_to_beat_best_cost() const noexcept { return cur_n_ - 2 * current_best_cost() + 2; }
    
    [[nodiscard]] bool gain_upper_bound_cannot_beat_best_cost(long double total_gain_ub) const noexcept {
        return as_int(std::floor(total_gain_ub + 1e-9L)) < gain_needed_to_beat_best_cost();
    }

    [[nodiscard]] std::span<const CompactPointIndex> line_points(int id) const noexcept {
        const std::size_t begin = line_point_offsets_[as_size(id)];
        const std::size_t end   = line_point_offsets_[as_size(id + 1)];
        return {line_point_pool_.data() + begin, end - begin};
    }

    [[nodiscard]] int positive_gain(const Task& task, int id) const noexcept {
        if (!task.available[as_size(id)]) return 0;
        const int gain = as_int(task.line_cover[as_size(id)]) - 2;
        return gain > 0 ? gain : 0;
    }

    [[nodiscard]] static UndoMark mark_undo(WorkerContext& ctx) {
        const std::size_t undo_token_depth_before = ctx.undo_token_stack.size();
        ctx.undo_token_stack.push_back(++ctx.next_undo_token);
        return { 
            ctx.disabled_lines.size(), 
            ctx.decremented_lines.size(), 
            ctx.newly_covered_points.size(),
            ctx.covered_point_degree.size(), 
            ctx.removed_from_productive.size(),
            ctx.removed_from_blocked_productive.size(),
            ctx.point_metadata_points.size(),
            undo_token_depth_before
        };
    }

    [[nodiscard]] std::size_t save_dual_for_points(WorkerContext& ctx, const std::vector<int>& points) const {
        const std::size_t mark = ctx.saved_dual_points.size();
        if (!ctx.lagrangian.dual_ready) return mark;
        for (const int p : points) {
            ctx.saved_dual_points.push_back(p);
            ctx.saved_dual_values.push_back(ctx.lagrangian.best_y[as_size(p)]);
        }
        return mark;
    }

    [[nodiscard]] std::size_t save_current_dual(WorkerContext& ctx) const { 
        return save_dual_for_points(ctx, ctx.lagrangian.active_points); 
    }

    void restore_saved_dual(WorkerContext& ctx, std::size_t mark) const {
        for (std::size_t i = mark; i < ctx.saved_dual_points.size(); ++i) {
            ctx.lagrangian.best_y[as_size(ctx.saved_dual_points[i])] = ctx.saved_dual_values[i];
        }
        ctx.saved_dual_points.resize(mark);
        ctx.saved_dual_values.resize(mark);
        ctx.lagrangian.dual_ready = true;
    }

    [[nodiscard]] static std::pair<int, int> root_productive_metrics(const Task& task) noexcept {
        int productive_incidence = 0;
        for (const int id : task.productive_ids) {
            productive_incidence += as_int(task.line_cover[as_size(id)]);
        }
        return { as_int(task.productive_ids.size()), productive_incidence };
    }

    void submit_candidate(int gain, const std::vector<int>& choice) {
        int observed = current_best_gain();
        while (gain > observed && !best_gain_.compare_exchange_weak(observed, gain, std::memory_order_relaxed, std::memory_order_relaxed)) {}
        
        if (gain < current_best_gain()) return;
        
        std::scoped_lock lock(best_choice_mutex_);
        if (gain > best_choice_gain_) {
            best_choice_gain_ = gain;
            best_choice_      = choice;
        }
    }

    // -----------------------------------------------------------------------
    // Core Private API Declarations
    // -----------------------------------------------------------------------
    [[nodiscard]] Task make_root_task(int n, const std::vector<int>& active_ids) const;
    void ensure_task_search_buffers(Task& task) const;
    void prepare_worker_context(WorkerContext& ctx) const;

    void dfs(Task& task, WorkerContext& ctx, int depth = 0);
    void run_task(Task task, WorkerContext& ctx);
    [[nodiscard]] std::vector<Task> build_frontier(Task root, WorkerContext& ctx);
    [[nodiscard]] SearchStats solve_frontier_parallel(const std::vector<Task>& frontier);

    [[nodiscard]] static bool better_branch_candidate(
        int gain_a, long double slack_a, long long pressure_a, int id_a,
        int gain_b, long double slack_b, long long pressure_b, int id_b) noexcept;
    [[nodiscard]] BranchChoice choose_branch(Task& task, WorkerContext& ctx, long double node_gain_ub) const;

    [[nodiscard]] long double lagrangian_upper_bound(const Task& task, WorkerContext& ctx, int max_iters) const;
    [[maybe_unused]] [[nodiscard]] long double lagrangian_upper_bound(const Task& task, WorkerContext& ctx) const;
    void restore_parent_lagrangian_view(const Task& task, WorkerContext& ctx) const;
    void build_compact_productive_lines(const Task& task, LagrangianScratch& scratch) const;
    void capture_root_cost_lb(long double gain_ub);

    [[nodiscard]] long double exclude_branch_upper_bound(const Task& task, const WorkerContext& ctx, int branch_line, long double node_gain_ub) const;
    [[nodiscard]] long double include_branch_upper_bound(Task& task, WorkerContext& ctx, int branch_line, long double node_gain_ub) const;

    template <typename SelectLineFn>
    void apply_forced_line_reduction(Task& task, WorkerContext& ctx, SelectLineFn&& select_line_fn, bool rebuild_metadata) const;
    void apply_dfs_forced_line_reduction(Task& task, WorkerContext& ctx) const;
    void apply_root_forced_line_reduction(Task& task, WorkerContext& ctx) const;

    void seed_from_line_ids(const Task& root, const std::vector<int>& line_ids);
    void greedy_complete_from_seed(const Task& base, WorkerContext& ctx, const std::vector<int>& line_ids);
    [[nodiscard]] int greedy_completion(const Task& task, WorkerContext& ctx) const;
    void maybe_submit_greedy(const Task& task, WorkerContext& ctx);

    void rebuild_productive_point_metadata(Task& task) const;
    void save_point_metadata_for_undo(const Task& task, WorkerContext& ctx, int p) const;
#if PC_VALIDATE_FORCE_METADATA
    [[nodiscard]] bool productive_point_metadata_matches_rebuild(const Task& task, int* mismatch_point = nullptr) const;
    void probe_validate_metadata(const Task& task, const char* stage) const;
#endif
    [[nodiscard]] int find_any_productive_line_for_point(const Task& task, int p) const;
    void on_productive_line_removed(Task& task, WorkerContext* ctx, int id) const;
    void on_productive_line_reinserted(Task& task, WorkerContext* ctx, int id) const;
    void remove_productive_impl(Task& task, WorkerContext* ctx, int id, std::vector<int>* removed_log) const;
    void remove_productive(Task& task, WorkerContext& ctx, int id) const;
    void reinsert_productive(Task& task, WorkerContext* ctx, int id) const;
    void add_blocked_productive_impl(Task& task, int id) const;
    void remove_blocked_productive_impl(Task& task, int id, std::vector<int>* removed_log) const;
    void remove_blocked_productive(Task& task, WorkerContext& ctx, int id) const;
    void disable_line_without_undo(Task& task, int id) const;
    void disable_line(Task& task, WorkerContext& ctx, int id) const;
    void cover_point(Task& task, WorkerContext& ctx, int p) const;
    void select_line(Task& task, WorkerContext& ctx, int id) const;
    void select_line_without_undo(Task& task, int id) const;
    void undo_to(Task& task, WorkerContext& ctx, UndoMark mark) const;

    [[nodiscard]] StateKey make_state_key(const Task& task, WorkerContext& ctx) const;
    [[nodiscard]] bool remember_frontier_state(const Task& task, WorkerContext& ctx, StateBestMap& best_gain) const;
    [[nodiscard]] bool frontier_state_is_stale(const Task& task, WorkerContext& ctx, const StateBestMap& best_gain) const;

    // -----------------------------------------------------------------------
    // Member Variables
    // -----------------------------------------------------------------------
    const std::vector<HeavyLine>&           lines_;
    const std::vector<std::vector<int>>&    incidence_;
    unsigned                                worker_count_{1};
    int                                     cur_n_{0};
    int                                     words_{0};
    std::vector<std::vector<int>>           ordered_point_lines_;
    std::size_t                             total_incidence_{0};
    std::vector<std::uint32_t>              line_point_offsets_;
    std::vector<CompactPointIndex>          line_point_pool_;
    
    std::atomic<int>                        best_gain_{0};
    std::atomic<bool>                       abort_{false};
    
    int                                     best_choice_gain_{0};
    int                                     root_cost_lb_{0};
    int                                     raw_root_cost_lb_{0};
    int                                     root_cost_lb_floor_{0};
    bool                                    root_lb_captured_{false};
    // Fractional part of the root Lagrangian gain bound; near 1.0 means lb_raw is close to its next tighter integer.
    double                                  root_ub_frac_{0.0};
    // Reused across sweep steps and fresh frontier tasks so sibling subproblems
    // do not cold-start the dual bound from zero every time.
    mutable std::vector<DualValue>          initial_dual_seed_;
    mutable bool                            initial_dual_seed_ready_{false};
    mutable bool                            root_dual_seed_captured_{false};
    
    std::vector<int>                        best_choice_;
    std::mutex                              best_choice_mutex_;
};

// ===========================================================================
// ExactGainSolver Implementation: Lifecycle & Main Entry
// ===========================================================================

ExactGainSolver::ExactGainSolver(const std::vector<HeavyLine>& lines, const std::vector<std::vector<int>>& incidence)
    : lines_(lines)
    , incidence_(incidence)
    , worker_count_(std::max(1U, std::thread::hardware_concurrency()))
    , total_incidence_(std::ranges::fold_left(lines, std::size_t{0},
          [](std::size_t total, const HeavyLine& line) { return total + line.points.size(); }))
    , line_point_offsets_(lines.size() + 1)
{
    line_point_pool_.reserve(total_incidence_);
    for (std::size_t id = 0; id < lines_.size(); ++id) {
        line_point_offsets_[id] = static_cast<std::uint32_t>(line_point_pool_.size());
        for (const int p : lines_[id].points) {
            line_point_pool_.push_back(static_cast<CompactPointIndex>(p));
        }
    }
    line_point_offsets_[lines_.size()] = static_cast<std::uint32_t>(line_point_pool_.size());
}

GainSolveResult ExactGainSolver::solve_for_n(
    int n, const std::vector<int>& active_ids, const std::vector<int>& warm_start_line_ids,
    const std::vector<double>& warm_dual_seed,
    const std::vector<int>& greedy_seed_line_ids, int root_cost_lb_floor)
{
    using namespace std::chrono_literals;

    cur_n_              = n;
    words_              = words_for_n(n);
    root_cost_lb_       = 0;
    raw_root_cost_lb_   = 0;
    root_cost_lb_floor_ = root_cost_lb_floor;
    root_lb_captured_   = false;
    root_ub_frac_       = 0.0;
    best_choice_gain_   = 0;
    
    best_gain_.store(0, std::memory_order_relaxed);
    abort_.store(false, std::memory_order_relaxed);
    
    best_choice_.clear();
    best_choice_.reserve(lines_.size());

    initial_dual_seed_.clear();
    if (!warm_dual_seed.empty()) {
        const std::size_t seed_count = std::min<std::size_t>(warm_dual_seed.size(), as_size(n));
        initial_dual_seed_.assign(warm_dual_seed.begin(), warm_dual_seed.begin() + static_cast<std::ptrdiff_t>(seed_count));
        initial_dual_seed_.resize(as_size(n), DualValue{0.0});
        initial_dual_seed_ready_ = true;
    } else {
        initial_dual_seed_ready_ = false;
    }
    root_dual_seed_captured_ = false;

    ordered_point_lines_.assign(as_size(n), {});
    for (int p = 0; p < n; ++p) {
        auto& ordered = ordered_point_lines_[as_size(p)];
        ordered.reserve(incidence_[as_size(p)].size());
        
        for (const int id : incidence_[as_size(p)]) {
            if (lines_[as_size(id)].activate_at <= n) {
                ordered.push_back(id);
            }
        }
        
        std::ranges::stable_sort(ordered, [&](int lhs, int rhs) {
            const auto sa = line_points(lhs).size();
            const auto sb = line_points(rhs).size();
            if (sa != sb) return sa > sb;
            return lhs < rhs;
        });
    }

    Task root = make_root_task(n, active_ids);
    WorkerContext root_ctx;
    prepare_worker_context(root_ctx);
    root_ctx.stats = {};
    
    seed_from_line_ids(root, warm_start_line_ids);
    seed_from_line_ids(root, greedy_seed_line_ids);
    apply_root_forced_line_reduction(root, root_ctx);
    
    const auto [root_productive_lines, root_productive_incidence] = root_productive_metrics(root);
    const int greedy_gain = greedy_completion(root, root_ctx);
    
    if (greedy_gain > current_best_gain()) {
        submit_candidate(greedy_gain, root_ctx.greedy_lines);
    }
    
    greedy_complete_from_seed(root, root_ctx, warm_start_line_ids);

    // Watchdog thread: sets `abort_` upon reaching the per-N time limit.
    // The flag is consumed via a relaxed load in `lagrangian_upper_bound()` 
    // to avoid overhead on the DFS hot path.
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

    // Snapshot forced count after the initial root pass; frontier-build and DFS forcing accumulate after this point.
    const int forced_at_dispatch = root_ctx.stats.forced;
    // Capture parallel frontier size; default is 1 for the single-thread path.
    int frontier_size_captured = 1;
    SearchStats total_stats = root_ctx.stats;
    
    if (worker_count_ <= 1) {
        run_task(std::move(root), root_ctx);
        total_stats = root_ctx.stats;
    } else {
        std::vector<Task> frontier = build_frontier(std::move(root), root_ctx);
        frontier_size_captured = as_int(frontier.size()); // Capture before tasks are consumed by workers or fallback.
        total_stats = root_ctx.stats;
        
        if (frontier.size() == 1) {
            run_task(std::move(frontier.front()), root_ctx);
            total_stats = root_ctx.stats;
        } else if (!frontier.empty()) {
            total_stats.merge_from(solve_frontier_parallel(frontier));
        }
    }

    watchdog_done.store(true, std::memory_order_relaxed);
    if (watchdog_thread.joinable()) {
        watchdog_thread.join();
    }

    GainSolveResult out;
    out.best_gain                 = current_best_gain();
    out.root_cost_lb              = root_cost_lb_;
    out.raw_root_cost_lb          = raw_root_cost_lb_;
    out.root_productive_lines     = root_productive_lines;
    out.root_productive_incidence = root_productive_incidence;
    out.forced                    = total_stats.forced;
    out.forced_root               = forced_at_dispatch;
    out.lag_iters                 = total_stats.lag_iters;
    out.polish_sweeps             = total_stats.polish_sweeps;
    out.depth_max                 = total_stats.depth_max;
    out.lag_prune                 = total_stats.lag_prune;
    out.strong_branch             = total_stats.strong_branch;
    out.nodes                     = total_stats.nodes;
    out.frontier_size             = frontier_size_captured;
    out.root_ub_frac              = root_ub_frac_;
    out.chosen_line_ids           = std::move(best_choice_);
    out.dual_seed                 = initial_dual_seed_;
    out.timed_out                 = abort_.load(std::memory_order_relaxed);
    
    return out;
}

// ===========================================================================
// ExactGainSolver Implementation: Core Search & Orchestration
// ===========================================================================

// Core Branch-and-Bound recursive search.
void ExactGainSolver::dfs(Task& task, WorkerContext& ctx, int depth) {
    if (root_interval_closed()) return;
    ctx.stats.depth_max = std::max(ctx.stats.depth_max, depth);
    
    const UndoMark node_undo            = mark_undo(ctx);
    const int start_gain                = task.current_gain;
    const std::size_t start_choice_size = task.current_choice.size();
    
    auto rollback = [&]() {
        task.current_choice.resize(start_choice_size);
        task.current_gain = start_gain;
        undo_to(task, ctx, node_undo);
    };

    // Apply the Exclusive Dependency Rule using a fresh metadata rebuild so
    // that staleness from parent branches cannot cause incorrect forcing.
    // All forced selections are logged in ctx and will be undone by rollback().
    apply_dfs_forced_line_reduction(task, ctx);
    if (root_interval_closed()) { rollback(); return; }

    ++ctx.stats.nodes;
    if (root_interval_closed()) {
        rollback();
        return;
    }
    
    const int lag_max_iters = std::min(64, 8 + task.active_count / 8);
    const long double ub    = lagrangian_upper_bound(task, ctx, lag_max_iters);
    capture_root_cost_lb(ub);
    
    const bool lag_killed = gain_upper_bound_cannot_beat_best_cost(ub);
    if (root_interval_closed() || lag_killed) {
        if (lag_killed) ++ctx.stats.lag_prune;
        rollback();
        return;
    }

    const BranchChoice branch = choose_branch(task, ctx, ub);
    const int branch_line     = branch.line;
    
    if (branch_line < 0) {
        submit_candidate(task.current_gain, task.current_choice);
        rollback(); 
        return;
    }

    const long double include_ub = include_branch_upper_bound(task, ctx, branch_line, ub);
    const long double exclude_ub = branch.exclude_ub;
    const bool include_pruned    = gain_upper_bound_cannot_beat_best_cost(include_ub);
    const bool exclude_pruned    = gain_upper_bound_cannot_beat_best_cost(exclude_ub);
    
    if (include_pruned && exclude_pruned) { 
        rollback(); 
        return; 
    }

    const std::size_t dual_mark = save_current_dual(ctx);

    if (!include_pruned) {
        const UndoMark undo = mark_undo(ctx);
        const int delta     = positive_gain(task, branch_line);
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

    if (exclude_pruned) { 
        rollback(); 
        return; 
    }

    {
        const UndoMark undo = mark_undo(ctx);
        disable_line(task, ctx, branch_line);
        
        dfs(task, ctx, depth + 1);
        
        restore_saved_dual(ctx, dual_mark);
        undo_to(task, ctx, undo);
    }

    rollback();
}

void ExactGainSolver::run_task(Task task, WorkerContext& ctx) {
    ensure_task_search_buffers(task);
    prepare_worker_context(ctx);
    if (root_interval_closed()) return;
    
    apply_root_forced_line_reduction(task, ctx);
    dfs(task, ctx);
}

// Performs initial BFS/Greedy expansion to generate independent tasks for parallel workers.
std::vector<ExactGainSolver::Task> ExactGainSolver::build_frontier(Task root, WorkerContext& ctx) {
    std::vector<Task> frontier;
    frontier.reserve(std::max<std::size_t>(4, as_size(worker_count_) * 8U));
    frontier.push_back(std::move(root));

    std::vector<Task> children;
    // =========================================================================
    // Frontier‑multiplier ladder – target = worker_count × multiplier.
    // =========================================================================
    // The solver uses worker_count = hardware_concurrency() (the number of CPU threads).
    // The frontier size = worker_count × multiplier. Larger multipliers increase
    // parallelism but also memory usage.
    //
    // -------------------------------------------------------------------------
    // Measurements from a Google Cloud c4d‑highcpu‑8 machine (8 vCPUs, 15 GB RAM)
    // running the solver up to N=943:
    //
    //   65,536 tasks (multiplier 8192U) → RSS ≈ 11.6 GB → ~187 KB per task.
    //  131,072 tasks (multiplier 16384U) → extrapolated RSS ≈ 23.4 GB → unsafe.
    //
    //   The safe maximum multiplier on this machine is 8192U.
    //   Multiplier 16384U would need ~23.4 GB and exceeds the 15 GB available.
    //
    // -------------------------------------------------------------------------
    // Per‑task memory varies greatly between different computers. For example:
    //   - Google Cloud c4d‑highcpu‑8 (15 GB): ~184 KB per task
    //   - Some personal computers (e.g., i9 9900K with 15 GB RAM) may show
    //     ~125 KB per task at frontier 8192, making 8192U unsafe (needs ~43 GB).
    //   Your results will likely be different – always measure.
    //
    // -------------------------------------------------------------------------
    // WHY THE LADDER USES 8192U FOR THE HARDEST PROBLEMS (on c4d‑highcpu‑8):
    //   Think of the search space as a construction project.
    //   - A small house (easy N) needs only a few workers – multipliers 4U to 128U.
    //   - A large building (moderately hard N) needs hundreds of workers – multiplier 2048U.
    //   - A massive skyscraper (cost 128–132) needs thousands of well‑coordinated
    //     workers and prefabricated sections – but only if you have enough space
    //     (RAM) to host them. On the tested cloud machine, 8192U (65k tasks) was
    //     that "ace in the hole" – enough to saturate all cores without OOM.
    //   On your hardware, the ideal multiplier may be higher or lower;
    //   always use the diagnostic command to find the safe maximum.
    //
    // -------------------------------------------------------------------------
    // IMPORTANT TOOL TO FIND YOUR SAFE MULTIPLIER:
    // Always measure using the diagnostic command (the "Frontier analyzer")
    // described in the GitHub repository:
    //   https://github.com/jespergran98/prime-line-cover
    // It reads the running solver's RSS and the most recent frontier size,
    // then tells you which multipliers are safe on your hardware.
    // Do not guess – always use the measured value.
    //
    // -------------------------------------------------------------------------
    // General guidelines for safe multipliers:
    //   - Multiplier 8192U is safe on many machines with at least 15 GB RAM,
    //     but not all – always verify with the diagnostic command.
    //   - Multiplier 16384U may be safe if your measured per‑task memory is low
    //     enough that the estimated memory usage fits within your available RAM.
    //     Use the diagnostic command to check before enabling it.
    //   - Higher multipliers (32768U, 131072U) are rarely practical on typical
    //     hardware; enable them only after thorough testing with the diagnostic command.
    //
    // -------------------------------------------------------------------------
    const unsigned frontier_multiplier =
    //  current_best_cost() >= 143 ? 131072U  // rarely feasible (weeks to months estimate)
    //: current_best_cost() >= 138 ? 32768U   // rarely feasible (days to weeks estimate)
    //: current_best_cost() >= 133 ? 16384U   // only if your diagnostic command says it is safe (hours to days estimate)
        current_best_cost() >= 128 ?  8192U   // safe on wr run - c4d‑highcpu‑8 (15 GB) (330–5545s)
      : current_best_cost() >= 125 ?  2048U   // safe on personal i9 9900k (133–1253s)
      : current_best_cost() >= 121 ?   512U   // safe on most machines (25–515s)
      : current_best_cost() >= 113 ?   128U   // safe (3.5–45s)
      : current_best_cost() >= 93  ?    16U   // safe (sub‑second)
                                   :     4U;  // trivial (microsecond to decisecond)

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
        apply_root_forced_line_reduction(task, ctx);
        if (frontier_state_is_stale(task, ctx, frontier_best)) continue;

        ctx.lagrangian.dual_ready = false;
        const int lag_max_iters   = std::min(64, 8 + task.active_count / 8);
        const long double ub      = lagrangian_upper_bound(task, ctx, lag_max_iters);
        capture_root_cost_lb(ub);
        
        if (root_interval_closed()) { frontier.clear(); break; }
        if (gain_upper_bound_cannot_beat_best_cost(ub)) {
            if (frontier.empty()) break;
            continue;
        }

        const BranchChoice branch = choose_branch(task, ctx, ub);
        const int branch_line     = branch.line;
        if (branch_line < 0) {
            submit_candidate(task.current_gain, task.current_choice);
            if (frontier.empty()) break;
            continue;
        }

        children.clear();
        children.reserve(2);

        const long double include_ub = include_branch_upper_bound(task, ctx, branch_line, ub);
        if (!gain_upper_bound_cannot_beat_best_cost(include_ub)) {
            Task include = task;
            ensure_task_search_buffers(include);
            const int delta = positive_gain(include, branch_line);
            if (delta > 0) {
                include.current_gain += delta;
                include.current_choice.push_back(branch_line);
                select_line_without_undo(include, branch_line);
                children.push_back(std::move(include));
            }
        }

        if (!gain_upper_bound_cannot_beat_best_cost(branch.exclude_ub)) {
            Task exclude = std::move(task); // Can safely move the last use
            ensure_task_search_buffers(exclude);
            disable_line_without_undo(exclude, branch_line);
            children.push_back(std::move(exclude));
        }

        for (Task& child : children) {
            // Use incremental forced reduction (no full metadata rebuild)
            apply_forced_line_reduction(child, ctx,
                [&](int best_line) { select_line_without_undo(child, best_line); },
                false);
            if (!remember_frontier_state(child, ctx, frontier_best)) {
                frontier.push_back(std::move(child));
            }
        }
    }

    std::ranges::sort(frontier, [](const Task& a, const Task& b) {
        if (a.active_count != b.active_count) return a.active_count > b.active_count;
        if (a.productive_ids.size() != b.productive_ids.size()) return a.productive_ids.size() > b.productive_ids.size();
        return a.current_gain < b.current_gain;
    });

    return frontier;
}

// Distributes the frontier tasks across the thread pool and aggregates statistics.
ExactGainSolver::SearchStats ExactGainSolver::solve_frontier_parallel(const std::vector<Task>& frontier) {
    const unsigned threads = std::min<unsigned>(worker_count_, static_cast<unsigned>(frontier.size()));
    if (threads <= 1) {
        WorkerContext ctx;
        run_task(frontier.front(), ctx);
        return ctx.stats;
    }

    std::atomic<std::size_t> next{0};
    std::mutex               stats_mutex;
    SearchStats              total_stats;
    std::vector<std::thread> workers;
    workers.reserve(as_size(threads - 1));

    auto worker = [&]() {
        WorkerContext ctx;
        while (true) {
            const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
            if (idx >= frontier.size()) break;
            run_task(frontier[idx], ctx);
        }
        std::scoped_lock lock(stats_mutex);
        total_stats.merge_from(ctx.stats);
    };

    for (unsigned i = 1; i < threads; ++i) {
        workers.emplace_back(worker);
    }
    worker();
    
    for (std::thread& thread : workers) {
        thread.join();
    }
    return total_stats;
}

// ===========================================================================
// ExactGainSolver Implementation: Branching Logic
// ===========================================================================

bool ExactGainSolver::better_branch_candidate(
    int gain_a, long double slack_a, long long pressure_a, int id_a,
    int gain_b, long double slack_b, long long pressure_b, int id_b) noexcept
{
    return (gain_a > gain_b) ||
           (gain_a == gain_b && slack_a > slack_b + 1e-18L) ||
           (gain_a == gain_b && std::abs(slack_a - slack_b) <= 1e-18L && pressure_a > pressure_b) ||
           (gain_a == gain_b && std::abs(slack_a - slack_b) <= 1e-18L && pressure_a == pressure_b && id_a < id_b);
}

ExactGainSolver::BranchChoice ExactGainSolver::choose_branch(Task& task, WorkerContext& ctx, long double node_gain_ub) const {
    struct Candidate {
        int         id{-1};
        int         gain{-1};
        long double slack{-1.0L};
        long long   pressure{-1};
        long double exclude_ub{0.0L};
    };

    std::array<Candidate, kStrongBranchShortlist> shortlist{};
    int shortlist_size = 0;

    auto branch_pressure = [&](Candidate& candidate) -> long long {
        if (candidate.pressure >= 0) return candidate.pressure;
        long long pressure = 0;
        for (const CompactPointIndex point_cp : line_points(candidate.id)) {
            const int p = static_cast<int>(point_cp);
            if (p >= cur_n_) break;
            if (task.active.test(p)) {
                pressure += static_cast<long long>(ordered_point_lines_[as_size(p)].size());
            }
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
        if (pos < kStrongBranchShortlist) { 
            shortlist[as_size(pos)] = candidate; 
            ++shortlist_size; 
        } else if (!better_shortlist_candidate(candidate, shortlist.back())) {
            return;
        } else { 
            shortlist.back() = candidate; 
            pos = kStrongBranchShortlist - 1; 
        }

        while (pos > 0 && better_shortlist_candidate(shortlist[as_size(pos)], shortlist[as_size(pos - 1)])) {
            std::swap(shortlist[as_size(pos)], shortlist[as_size(pos - 1)]);
            --pos;
        }
    };

    const auto& line_slack = ctx.lagrangian.work_line_slack_compact;
    for (const int id : task.productive_ids) {
        const int li = task.productive_pos[as_size(id)];
        insert_candidate({
            id,
            as_int(task.line_cover[as_size(id)]) - 2,
            std::max(0.0L, static_cast<long double>(line_slack[as_size(li)])),
            -1,
            0.0L
        });
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
    ++ctx.stats.strong_branch;

    const std::vector<int> parent_active_points = ctx.lagrangian.active_points;
    int best_pruned_children     = -1;
    long double best_worst_child_ub = std::numeric_limits<long double>::infinity();
    long double best_total_child_ub = std::numeric_limits<long double>::infinity();
    int best_gain                = -1;
    long double best_slack          = -1.0L;
    long long best_pressure      = -1;

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
            include_ub = lagrangian_upper_bound(task, ctx, kStrongBranchMaxIters);
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
            best_pruned_children = pruned_children; 
            best_worst_child_ub  = worst_child_ub; 
            best_total_child_ub  = total_child_ub;
            best_gain            = candidate.gain; 
            best_slack           = candidate.slack; 
            best_pressure        = candidate_pressure;
            best.line            = candidate.id; 
            best.exclude_ub      = candidate.exclude_ub;
        }
    }
    
    restore_parent_lagrangian_view(task, ctx);
    return best;
}

// ===========================================================================
// ExactGainSolver Implementation: Bounding & Lagrangian
// ===========================================================================

long double ExactGainSolver::lagrangian_upper_bound(const Task& task, WorkerContext& ctx, int max_iters) const {
    // Fast cancellation check. Returning -1 forces immediate node pruning in DFS.
    if (abort_.load(std::memory_order_relaxed)) return -1.0L;

    auto& scratch            = ctx.lagrangian;
    auto& active_y           = scratch.active_y;
    auto& active_best_y      = scratch.active_best_y;
    auto& best_y             = scratch.best_y;
    auto& active_subgrad     = scratch.active_subgrad;
    auto& compact_line_slack = scratch.work_line_slack_compact;

    build_compact_productive_lines(task, scratch);
    
    const auto& active_points            = scratch.active_points;
    const auto& compact_line_ids         = scratch.compact_line_ids;
    const auto& compact_offsets          = scratch.compact_line_offsets;
    const auto& compact_point_aps        = scratch.compact_point_aps;
    const auto& compact_eval3            = scratch.compact_eval3;
    const auto& compact_eval4            = scratch.compact_eval4;
    const auto& compact_eval5            = scratch.compact_eval5;
    const auto& compact_eval_other_lines = scratch.compact_eval_other_lines;
    const auto& compact_rem2             = scratch.compact_rem2;
    
    const int active_point_count         = as_int(active_points.size());
    const auto* compact_offsets_data     = compact_offsets.data();
    const auto* compact_point_aps_data   = compact_point_aps.data();
    const auto* compact_rem2_data        = compact_rem2.data();
    const bool warm_start_dual           = scratch.dual_ready;

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
        std::fill_n(active_y.data(), active_point_count, DualValue{0.0});
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            const int rem2 = compact_rem2_data[li];
            const int rem  = rem2 + 2;
            const DualValue density = static_cast<DualValue>(rem2) / static_cast<DualValue>(rem);
            
            for (std::size_t pos = compact_offsets_data[li]; pos < compact_offsets_data[li + 1]; ++pos) {
                DualValue& slot = active_y[as_size(compact_point_aps_data[pos])];
                if (density > slot) slot = density;
            }
        }
        for (int ai = 0; ai < active_point_count; ++ai) {
            active_y[as_size(ai)] = std::min<DualValue>(1.0, active_y[as_size(ai)]);
        }
    }

    auto eval_dense_fast = [&](const std::vector<DualValue>& yy, std::vector<SubgradientValue>* subgradient) -> DualValue {
        const auto* yy_data      = yy.data();
        auto* subgradient_data   = subgradient ? subgradient->data() : nullptr;
        DualValue value          = 0.0;
        
        for (int ai = 0; ai < active_point_count; ++ai) {
            value += yy_data[ai];
            if (subgradient_data) subgradient_data[ai] = SubgradientValue{1};
        }
        
        for (const CompactLine3& line : compact_eval3) {
            const DualValue rc = 1.0 - (yy_data[line.a] + yy_data[line.b] + yy_data[line.c]);
            if (rc <= 1e-18) continue;
            value += rc;
            if (subgradient_data) {
                --subgradient_data[line.a];
                --subgradient_data[line.b];
                --subgradient_data[line.c];
            }
        }
        
        for (const CompactLine4& line : compact_eval4) {
            const DualValue rc = 2.0 - (yy_data[line.a] + yy_data[line.b] + yy_data[line.c] + yy_data[line.d]);
            if (rc <= 1e-18) continue;
            value += rc;
            if (subgradient_data) {
                --subgradient_data[line.a];
                --subgradient_data[line.b];
                --subgradient_data[line.c];
                --subgradient_data[line.d];
            }
        }
        
        for (const CompactLine5& line : compact_eval5) {
            const DualValue rc = 3.0 - (yy_data[line.a] + yy_data[line.b] + yy_data[line.c] + yy_data[line.d] + yy_data[line.e]);
            if (rc <= 1e-18) continue;
            value += rc;
            if (subgradient_data) {
                --subgradient_data[line.a];
                --subgradient_data[line.b];
                --subgradient_data[line.c];
                --subgradient_data[line.d];
                --subgradient_data[line.e];
            }
        }
        
        for (const int li_int : compact_eval_other_lines) {
            const std::size_t li    = as_size(li_int);
            const std::size_t begin = compact_offsets_data[li];
            const std::size_t end   = compact_offsets_data[li + 1];
            DualValue sum_y         = 0.0;
            
            for (std::size_t pos = begin; pos < end; ++pos) {
                sum_y += yy_data[compact_point_aps_data[pos]];
            }
                
            const DualValue rc = static_cast<DualValue>(compact_rem2_data[li]) - sum_y;
            if (rc > 1e-18) {
                value += rc;
                if (subgradient_data) {
                    for (std::size_t pos = begin; pos < end; ++pos) {
                        --subgradient_data[compact_point_aps_data[pos]];
                    }
                }
            }
        }
        return value;
    };

    auto eval_dense_exact = [&](const std::vector<DualValue>& yy, bool store_slack) -> long double {
        const auto* yy_data = yy.data();
        long double value   = 0.0L;
        
        for (int ai = 0; ai < active_point_count; ++ai) {
            value += static_cast<long double>(yy_data[ai]);
        }
            
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            const std::size_t begin = compact_offsets_data[li];
            const std::size_t end   = compact_offsets_data[li + 1];
            long double sum_y       = 0.0L;
            
            for (std::size_t pos = begin; pos < end; ++pos) {
                sum_y += static_cast<long double>(yy_data[compact_point_aps_data[pos]]);
            }
                
            const long double rc = static_cast<long double>(compact_rem2_data[li]) - sum_y;
            if (store_slack) compact_line_slack[li] = static_cast<DualValue>(rc);
            if (rc > 1e-18L) value += rc;
        }
        return value;
    };

    auto store_best_dual = [&]() {
        for (int ai = 0; ai < active_point_count; ++ai) {
            best_y[as_size(active_points[as_size(ai)])] = active_best_y[as_size(ai)];
        }
        if (!root_dual_seed_captured_) {
            initial_dual_seed_.assign(as_size(cur_n_), DualValue{0.0});
            for (int ai = 0; ai < active_point_count; ++ai) {
                initial_dual_seed_[as_size(active_points[as_size(ai)])] = active_best_y[as_size(ai)];
            }
            initial_dual_seed_ready_ = true;
            root_dual_seed_captured_ = true;
        }
    };

    const long double best_exact = eval_dense_exact(active_y, false);
    DualValue best               = static_cast<DualValue>(best_exact);
    
    for (int ai = 0; ai < active_point_count; ++ai) {
        active_best_y[as_size(ai)] = active_y[as_size(ai)];
    }

    const long double initial_total = static_cast<long double>(task.current_gain) + best_exact;
    scratch.dual_ready = true;
    
    if (gain_upper_bound_cannot_beat_best_cost(initial_total)) {
        store_best_dual();
        return initial_total;
    }

    long double confirmed_prune_total = 0.0L;
    const bool aggressive_tail_warm = warm_start_dual &&
                                      current_best_cost() >= 113 &&
                                      task.current_choice.size() >= 4;
    auto maybe_confirm_prune = [&](DualValue candidate_best) -> bool {
        const long double optimistic_total = static_cast<long double>(task.current_gain) + static_cast<long double>(candidate_best);
        if (optimistic_total + 1e-9L >= static_cast<long double>(gain_needed_to_beat_best_cost())) return false;
        
        const long double exact_total = static_cast<long double>(task.current_gain) + eval_dense_exact(active_best_y, false);
        if (!gain_upper_bound_cannot_beat_best_cost(exact_total)) return false;
        
        confirmed_prune_total = exact_total;
        store_best_dual();
        return true;
    };

    auto run_subgradient_iters = [&](int begin_iter, int end_iter) -> bool {
        for (int iter = begin_iter; iter < end_iter; ++iter) {
            ++ctx.stats.lag_iters;
            const DualValue value = eval_dense_fast(active_y, &active_subgrad);
            
            if (value + 1e-18 < best) {
                best = value;
                for (int ai = 0; ai < active_point_count; ++ai) {
                    active_best_y[as_size(ai)] = active_y[as_size(ai)];
                }
                if (maybe_confirm_prune(best)) return true;
            }

            DualValue norm2 = 0.0;
            for (int ai = 0; ai < active_point_count; ++ai) {
                const DualValue g = static_cast<DualValue>(active_subgrad[as_size(ai)]);
                norm2 += g * g;
            }
            if (norm2 <= 1e-20) return false;

            const DualValue incumbent_gap = static_cast<DualValue>(std::max(0, gain_needed_to_beat_best_cost() - task.current_gain));
            const DualValue target        = incumbent_gap;
            const DualValue gap           = std::max<DualValue>(0.0, value - target);
            if (gap <= 1e-9) return false;

            const DualValue step_scale = aggressive_tail_warm ? DualValue{1.75}
                                       : warm_start_dual      ? DualValue{1.5}
                                                              : DualValue{1.35};
            const DualValue step = (step_scale * gap) / norm2;
            for (int ai = 0; ai < active_point_count; ++ai) {
                const DualValue g = static_cast<DualValue>(active_subgrad[as_size(ai)]);
                DualValue next    = active_y[as_size(ai)] - step * g;
                
                if (next < 0.0) next = 0.0;
                else if (next > 1.0) next = 1.0;
                
                active_y[as_size(ai)] = next;
            }
        }
        return false;
    };

    constexpr int kWarmStartFirstStageIters    = 16;
    // Warm-started duals already capture most of the useful structure on the
    // gap-1 tail, so only pay for the full second-stage refinement when the
    // 16-iteration pass is genuinely close to a pruning threshold.
    constexpr long double kWarmStartRefineWindow = 1.0L;
    const bool use_short_warm_stage = warm_start_dual &&
                                      current_best_cost() >= 112 &&
                                      task.current_choice.size() >= 4;
                                      
    const int first_stage_iters = use_short_warm_stage ? std::min(max_iters, kWarmStartFirstStageIters) : max_iters;
    if (run_subgradient_iters(0, first_stage_iters)) return confirmed_prune_total;
    
    if (use_short_warm_stage && first_stage_iters < max_iters) {
        const long double staged_total = static_cast<long double>(task.current_gain) + static_cast<long double>(best);
        if (staged_total < static_cast<long double>(gain_needed_to_beat_best_cost()) + kWarmStartRefineWindow) {
            if (run_subgradient_iters(first_stage_iters, max_iters)) return confirmed_prune_total;
        }
    }

    if (maybe_confirm_prune(best)) return confirmed_prune_total;

    const int kCoordinatePolishMaxSweeps = scratch.dual_ready ? 1 : 3;
    constexpr long double kCoordinatePolishWindow = 0.5L;
    const long double best_total = static_cast<long double>(task.current_gain) + static_cast<long double>(best);
    
    // Gate coordinate polish to nodes near the pruning threshold, as the 
    // single-sweep polish consumes a notable portion of bound time.
    if (kCoordinatePolishMaxSweeps > 0 &&
        best_total < static_cast<long double>(gain_needed_to_beat_best_cost()) + kCoordinatePolishWindow) 
    {
        const int ap_size = active_point_count;
        auto& ap_offsets  = scratch.compact_pt_prod_offsets;
        auto& ap_flat     = scratch.compact_pt_prod_flat;
        auto& fill_pos    = scratch.compact_ap_fill_pos;

        ap_offsets.assign(as_size(ap_size + 1), 0);
        for (const CompactPointIndex ai_cp : compact_point_aps) {
            ++ap_offsets[as_size(static_cast<int>(ai_cp) + 1)];
        }

        for (int i = 0; i < ap_size; ++i) {
            ap_offsets[as_size(i + 1)] += ap_offsets[as_size(i)];
        }

        ap_flat.resize(as_size(ap_offsets[as_size(ap_size)]));
        fill_pos.resize(as_size(ap_size));
        for (int i = 0; i < ap_size; ++i) {
            fill_pos[as_size(i)] = ap_offsets[as_size(i)];
        }

        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            for (std::size_t pos = compact_offsets_data[li]; pos < compact_offsets_data[li + 1]; ++pos) {
                const int ai = compact_point_aps_data[pos];
                ap_flat[as_size(fill_pos[as_size(ai)]++)] = as_int(li);
            }
        }

        auto& wls_c = scratch.work_line_slack_compact;
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            const std::size_t begin = compact_offsets_data[li];
            const std::size_t end   = compact_offsets_data[li + 1];
            DualValue slack         = static_cast<DualValue>(compact_rem2_data[li]);
            
            for (std::size_t pos = begin; pos < end; ++pos) {
                slack -= active_best_y[compact_point_aps_data[pos]];
            }
            wls_c[li] = slack;
        }

        auto slack_objective = [&]() -> DualValue {
            DualValue value = 0.0;
            for (int ai = 0; ai < active_point_count; ++ai) {
                value += active_best_y[as_size(ai)];
            }
            for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
                if (wls_c[li] > 1e-18) value += wls_c[li];
            }
            return value;
        };

        DualValue polished_best = best;
        for (int sweep = 0; sweep < kCoordinatePolishMaxSweeps; ++sweep) {
            ++ctx.stats.polish_sweeps; // Count each polish sweep for bound-kernel diagnostics.
            bool changed = false;
            for (int ai = 0; ai < active_point_count; ++ai) {
                const DualValue old = active_best_y[as_size(ai)];
                DualValue largest   = 0.0;
                DualValue second    = 0.0;
                const int beg       = ap_offsets[as_size(ai)];
                const int end       = ap_offsets[as_size(ai + 1)];
                
                for (int k = beg; k < end; ++k) {
                    const int li = ap_flat[as_size(k)];
                    const DualValue adjusted = wls_c[as_size(li)] + old;
                    
                    if (adjusted <= 1e-18) continue;
                    if (adjusted >= largest) { 
                        second = largest; 
                        largest = adjusted; 
                    } else if (adjusted > second) { 
                        second = adjusted; 
                    }
                }

                const DualValue next = std::min<DualValue>(1.0, second);
                if (std::abs(next - old) <= 1e-15) continue;
                
                changed = true;
                active_best_y[as_size(ai)] = next;
                const DualValue delta = next - old;
                for (int k = beg; k < end; ++k) {
                    wls_c[as_size(ap_flat[as_size(k)])] -= delta;
                }
            }
            polished_best = std::min(polished_best, slack_objective());
            if (!changed) break;
        }
        best = polished_best;
    }

    store_best_dual();
    return static_cast<long double>(task.current_gain) + eval_dense_exact(active_best_y, true);
}

[[maybe_unused]] long double ExactGainSolver::lagrangian_upper_bound(const Task& task, WorkerContext& ctx) const {
    return lagrangian_upper_bound(task, ctx, 64); 
}

void ExactGainSolver::restore_parent_lagrangian_view(const Task& task, WorkerContext& ctx) const {
    auto& scratch = ctx.lagrangian;
    if (!scratch.dual_ready) return;

    build_compact_productive_lines(task, scratch);
    
    const auto& active_points     = scratch.active_points;
    const auto& compact_offsets   = scratch.compact_line_offsets;
    const auto& compact_point_aps = scratch.compact_point_aps;
    const auto& compact_rem2      = scratch.compact_rem2;
    const auto& best_y            = scratch.best_y;
    auto& line_slack              = scratch.work_line_slack_compact;

    for (std::size_t li = 0; li < scratch.compact_line_ids.size(); ++li) {
        const std::size_t begin = compact_offsets[li];
        const std::size_t end   = compact_offsets[li + 1];
        DualValue slack         = static_cast<DualValue>(compact_rem2[li]);
        
        for (std::size_t pos = begin; pos < end; ++pos) {
            const int p = active_points[as_size(compact_point_aps[pos])];
            slack -= best_y[as_size(p)];
        }
        line_slack[li] = slack;
    }
}

void ExactGainSolver::build_compact_productive_lines(const Task& task, LagrangianScratch& scratch) const {
    auto& active_points            = scratch.active_points;
    auto& compact_line_ids         = scratch.compact_line_ids;
    auto& compact_offsets          = scratch.compact_line_offsets;
    auto& compact_points           = scratch.compact_line_points;
    auto& compact_point_aps        = scratch.compact_point_aps;
    auto& compact_eval3            = scratch.compact_eval3;
    auto& compact_eval4            = scratch.compact_eval4;
    auto& compact_eval5            = scratch.compact_eval5;
    auto& compact_eval_other_lines = scratch.compact_eval_other_lines;

    active_points.clear();
    compact_line_ids.clear(); 
    compact_offsets.clear(); 
    compact_points.clear(); 
    compact_point_aps.clear();
    compact_eval3.clear(); 
    compact_eval4.clear(); 
    compact_eval5.clear(); 
    compact_eval_other_lines.clear();
    
    const std::size_t productive_count = task.productive_ids.size();
    
    if (compact_line_ids.capacity() < productive_count) compact_line_ids.reserve(productive_count);
    if (compact_offsets.capacity() < productive_count + 1) compact_offsets.reserve(productive_count + 1);

    std::size_t total_active_incidence = 0;
    for (const int id : task.productive_ids) {
        total_active_incidence += as_size(task.line_cover[as_size(id)]);
    }
    if (compact_points.capacity() < total_active_incidence) compact_points.reserve(total_active_incidence);
    if (compact_point_aps.capacity() < total_active_incidence) compact_point_aps.reserve(total_active_incidence);

    BitMask1024 productive_union;
    compact_offsets.push_back(0);
    for (const int id : task.productive_ids) {
        compact_line_ids.push_back(id);
        int remaining = as_int(task.line_cover[as_size(id)]);
        
        for (const CompactPointIndex point_cp : line_points(id)) {
            const int p = static_cast<int>(point_cp);
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            
            productive_union.set(p);
            compact_points.push_back(static_cast<CompactPointIndex>(p));
            if (--remaining == 0) break;
        }
        compact_offsets.push_back(compact_points.size());
    }

    // Recover the active-point set simultaneously to avoid a redundant global scan.
    productive_union.for_each_set_bit(words_, [&](int p) { active_points.push_back(p); });

    {
        // Precompute `line_cover[id] - 2` to replace random access in `eval_dense` loops 
        // with sequential reads. This drastically reduces the L1 cache footprint.
        auto& rem2 = scratch.compact_rem2;
        rem2.resize(compact_line_ids.size());
        for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
            rem2[li] = as_int(task.line_cover[as_size(compact_line_ids[li])]) - 2;
        }
    }

    const int ap_size = as_int(active_points.size());
    auto& pt_to_ap    = scratch.point_to_ap_scratch;

    for (int ai = 0; ai < ap_size; ++ai) pt_to_ap[as_size(active_points[as_size(ai)])] = ai;

    compact_point_aps.resize(compact_points.size());
    for (std::size_t pos = 0; pos < compact_points.size(); ++pos) {
        compact_point_aps[pos] = static_cast<CompactPointIndex>(pt_to_ap[as_size(compact_points[pos])]);
    }

    for (const int p : active_points) pt_to_ap[as_size(p)] = -1;

    compact_eval3.reserve(compact_line_ids.size());
    compact_eval4.reserve(compact_line_ids.size());
    compact_eval5.reserve(compact_line_ids.size());
    compact_eval_other_lines.reserve(compact_line_ids.size());
    
    for (std::size_t li = 0; li < compact_line_ids.size(); ++li) {
        const std::size_t begin = compact_offsets[li];
        const std::size_t end   = compact_offsets[li + 1];
        const std::size_t rem   = end - begin;
        
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

    scratch.work_line_slack_compact.resize(compact_line_ids.size());
}

void ExactGainSolver::capture_root_cost_lb(long double gain_ub) {
    if (root_lb_captured_) return;
    raw_root_cost_lb_ = cost_from_gain(cur_n_, as_int(std::floor(gain_ub + 1e-9L)));
    root_cost_lb_     = std::max(raw_root_cost_lb_, root_cost_lb_floor_);
    root_lb_captured_ = true;
    // Diagnostic: fractional part of the gain bound, independent of the algorithmic floor above.
    const double gain_ub_d = static_cast<double>(gain_ub);
    root_ub_frac_ = gain_ub_d - std::floor(gain_ub_d + 1e-9);
}

long double ExactGainSolver::exclude_branch_upper_bound(const Task& task, const WorkerContext& ctx, int branch_line, long double node_gain_ub) const {
    const auto& best_y     = ctx.lagrangian.best_y;
    const auto& line_slack = ctx.lagrangian.work_line_slack_compact;
    const int branch_li    = task.productive_pos[as_size(branch_line)];
    long double capped     = node_gain_ub - std::max(0.0L, static_cast<long double>(line_slack[as_size(branch_li)]));

    for (const CompactPointIndex point_cp : line_points(branch_line)) {
        const int p = static_cast<int>(point_cp);
        if (p >= cur_n_) break;
        if (!task.active.test(p)) continue;
        
        const long double residual = 1.0L - best_y[as_size(p)];
        if (residual <= 1e-18L) continue;

        long double requiring_sum  = 0.0L;
        long double requiring_best = 0.0L;
        
        for (const int id : ordered_point_lines_[as_size(p)]) {
            if (id == branch_line || !task.available[as_size(id)] || task.line_cover[as_size(id)] < 3) continue;
            
            const int li         = task.productive_pos[as_size(id)];
            const long double rc = static_cast<long double>(line_slack[as_size(li)]);
            
            if (rc <= 1e-18L || rc > residual + 1e-18L) continue;
            requiring_sum += rc;
            requiring_best = std::max(requiring_best, rc);
        }
        capped -= (requiring_sum - requiring_best);
    }
    return std::max<long double>(static_cast<long double>(task.current_gain), capped);
}

long double ExactGainSolver::include_branch_upper_bound(Task& task, WorkerContext& ctx, int branch_line, long double node_gain_ub) const {
    const auto& best_y     = ctx.lagrangian.best_y;
    const auto& line_slack = ctx.lagrangian.work_line_slack_compact;
    auto& touched          = ctx.include_ub_touched_lines;
    auto& reduction        = ctx.include_ub_line_reduction;
    touched.clear();

    const int delta     = positive_gain(task, branch_line);
    const int branch_li = task.productive_pos[as_size(branch_line)];
    long double capped  = node_gain_ub + std::min<long double>(0.0L, static_cast<long double>(line_slack[as_size(branch_li)]));

    for (const CompactPointIndex point_cp : line_points(branch_line)) {
        const int p = static_cast<int>(point_cp);
        if (p >= cur_n_) break;
        if (!task.active.test(p)) continue;
        
        const DualValue residual = 1.0 - best_y[as_size(p)];
        if (residual <= 1e-18) continue;

        for (const int id : ordered_point_lines_[as_size(p)]) {
            if (id == branch_line || !task.available[as_size(id)] || task.line_cover[as_size(id)] < 3) continue;
            
            const int li = task.productive_pos[as_size(id)];
            if (li < 0) continue;
            
            if (reduction[as_size(li)] == 0.0) touched.push_back(li);
            reduction[as_size(li)] += residual;
        }
    }

    for (const int li : touched) {
        const long double old_rc = static_cast<long double>(line_slack[as_size(li)]);
        const long double new_rc = old_rc - static_cast<long double>(reduction[as_size(li)]);
        
        capped += std::max<long double>(0.0L, new_rc) - std::max<long double>(0.0L, old_rc);
        reduction[as_size(li)] = 0.0;
    }

    return std::max<long double>(static_cast<long double>(task.current_gain + delta), capped);
}

// ===========================================================================
// ExactGainSolver Implementation: Heuristics & Reductions
// ===========================================================================

//
// Exclusive Dependency Rule (proved, unconditional)
//
// Let L be a productive heavy line (covers at least 3 active points).
// Let U be the set of active points for which L is their *only* available
// productive heavy line (no other heavy line with coverage >= 3 contains
// any point of U). The rule applies when |U| >= 3. The proof shows delta <= 0
// for |U| = 3 (delta = 0) and delta < 0 for |U| >= 4, so forcing L never worsens
// the solution.
//
// Proof that L belongs to some optimal completion:
//   Consider any optimal solution that omits L.  The points of U are
//   covered by lines that are not productive heavy.  Each such line
//   covers at most 2 active points in total, and at most 1 of them can
//   lie outside U (because it must cover at least one point of U).
//   Let k be the number of those lines, and let d <= k be the number of
//   non-U points they cover.  Replace these k lines with L (cost 1) and
//   cover the d non-U points using at most ceil(d/2) <= ceil(k/2) new pairing
//   lines.  The change in total cost is
//        delta = -k + 1 + ceil(k/2).
//   For k = 2, delta = 0; for k = 3, delta = 0; for k >= 4, delta <= -k/2 + 1.5 < 0.
//   Because |U| >= 3 forces k >= ceil(|U|/2) >= 2, we have delta <= 0.  Hence there
//   exists an optimal solution that includes L, and forcing L never
//   worsens the solution.  The rule is unconditional.
//
// Consequently L may be forced into the partial solution without any
// loss of optimality.  The threshold is 3.
//
template <typename SelectLineFn>
void ExactGainSolver::apply_forced_line_reduction(Task& task, WorkerContext& ctx, SelectLineFn&& select_line_fn, bool rebuild_metadata) const {
    auto& touched = ctx.forced_touched_lines;

    while (true) {
#if PC_VALIDATE_FORCE_METADATA
        probe_validate_metadata(task, "apply_forced_line_reduction");
#endif
#if !PC_SKIP_FORCE_REBUILD_PRODUCTIVE_METADATA
        rebuild_productive_point_metadata(task);
#else
        if (rebuild_metadata) rebuild_productive_point_metadata(task);
#endif

        touched.clear();
        task.active.for_each_set_bit(words_, [&](int p) {
            const int unique_line = task.sole_productive_line[as_size(p)];
            if (unique_line >= 0 && (touched.empty() || touched.back() != unique_line)) {
                touched.push_back(unique_line);
            }
        });

        int best_line = -1;
        int best_unique = 2;
        int best_gain = -1;
        for (const int id : touched) {
            int unique_points = 0;
            for (const CompactPointIndex cp : line_points(id)) {
                const int p = static_cast<int>(cp);
                if (p >= cur_n_) break;
                if (!task.active.test(p)) continue;
                if (task.sole_productive_line[as_size(p)] == id) ++unique_points;
            }
            if (unique_points < 3) continue;

            const int gain = positive_gain(task, id);
            if (gain <= 0) continue;

            const bool better =
                (unique_points > best_unique) ||
                (unique_points == best_unique && gain > best_gain) ||
                (unique_points == best_unique && gain == best_gain &&
                 (best_line < 0 || id < best_line));
            if (better) {
                best_line   = id;
                best_unique = unique_points;
                best_gain   = gain;
            }
        }
        if (best_line < 0) break;

        ++ctx.stats.forced;
        select_line_fn(best_line);
        task.current_gain += best_gain;
        task.current_choice.push_back(best_line);
    }
}

// DFS-safe variant: uses the undo-aware selector so forced choices roll
// back cleanly with the surrounding DFS node.
void ExactGainSolver::apply_dfs_forced_line_reduction(Task& task, WorkerContext& ctx) const {
    apply_forced_line_reduction(task, ctx, [&](int best_line) {
        select_line(task, ctx, best_line);
    }, false);
}

// Root/frontier variant: uses the non-undo selector because the caller
// owns the task outright and no DFS rollback is required.
void ExactGainSolver::apply_root_forced_line_reduction(Task& task, WorkerContext& ctx) const {
    apply_forced_line_reduction(task, ctx, [&](int best_line) {
        select_line_without_undo(task, best_line);
    }, true);
}

int ExactGainSolver::greedy_completion(const Task& task, WorkerContext& ctx) const {
    BitMask1024 temp_active = task.active;
    const int max_gain      = std::max(0, task.active_count - 2);
    
    std::fill_n(ctx.greedy_bucket_head.begin(), max_gain + 1, -1);
    
    for (const int id : task.productive_ids) {
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
            for (const CompactPointIndex point_cp : line_points(id)) {
                const int p = static_cast<int>(point_cp);
                if (p >= cur_n_) break;
                if (temp_active.test(p)) ++current_gain;
            }
            if (current_gain <= 0) continue;

            gain += current_gain;
            chosen.push_back(id);
            
            for (const CompactPointIndex point_cp : line_points(id)) {
                const int p = static_cast<int>(point_cp);
                if (p >= cur_n_) break;
                if (temp_active.test(p)) temp_active.reset(p);
            }
        }
    }
    return gain;
}

void ExactGainSolver::seed_from_line_ids(const Task& root, const std::vector<int>& line_ids) {
    if (line_ids.empty()) return;
    
    Task trial = root;
    int gain = 0;
    std::vector<int> chosen;
    chosen.reserve(line_ids.size());
    
    for (const int id : line_ids) {
        if (id < 0 || id >= as_int(lines_.size()) || lines_[as_size(id)].activate_at > cur_n_) continue;
        const int delta = positive_gain(trial, id);
        if (delta <= 0) continue;
        
        gain += delta;
        chosen.push_back(id);
        select_line_without_undo(trial, id);
    }
    if (gain > current_best_gain()) submit_candidate(gain, chosen);
}

void ExactGainSolver::greedy_complete_from_seed(const Task& base, WorkerContext& ctx, const std::vector<int>& line_ids) {
    if (line_ids.empty()) return;
    
    Task trial = base;
    ensure_task_search_buffers(trial);
    
    for (const int id : line_ids) {
        if (id < 0 || id >= as_int(lines_.size()) || lines_[as_size(id)].activate_at > cur_n_) continue;
        const int delta = positive_gain(trial, id);
        if (delta <= 0) continue;
        
        select_line_without_undo(trial, id);
        trial.current_gain += delta;
        trial.current_choice.push_back(id);
    }
    maybe_submit_greedy(trial, ctx);
}

void ExactGainSolver::maybe_submit_greedy(const Task& task, WorkerContext& ctx) {
    const int greedy_gain = greedy_completion(task, ctx);
    if (task.current_gain + greedy_gain <= current_best_gain()) return;

    ctx.candidate_choice = task.current_choice;
    ctx.candidate_choice.insert(ctx.candidate_choice.end(), ctx.greedy_lines.begin(), ctx.greedy_lines.end());
    submit_candidate(task.current_gain + greedy_gain, ctx.candidate_choice);
}

// ===========================================================================
// ExactGainSolver Implementation: Task State Management
// ===========================================================================

ExactGainSolver::Task ExactGainSolver::make_root_task(int n, const std::vector<int>& active_ids) const {
    Task root;
    root.active       = BitMask1024::prefix(n);
    root.active_count = n;
    
    root.line_cover.assign(lines_.size(), 0);
    root.productive_degree.assign(as_size(n), 0);
    root.sole_productive_line.assign(as_size(n), kNoProductiveLine);
    root.available.assign(lines_.size(), 0);
    root.productive_pos.assign(lines_.size(), -1);
    root.blocked_productive_pos.assign(lines_.size(), -1);
    
    root.current_choice.reserve(lines_.size());
    root.productive_ids.reserve(active_ids.size());
    root.blocked_productive_ids.reserve(64);
    
    for (const int id : active_ids) {
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

void ExactGainSolver::ensure_task_search_buffers(Task& task) const {
    if (task.current_choice.capacity() < lines_.size()) task.current_choice.reserve(lines_.size());
    if (task.blocked_productive_ids.capacity() < 64) task.blocked_productive_ids.reserve(64);
}

void ExactGainSolver::prepare_worker_context(WorkerContext& ctx) const {
    ctx.lagrangian.ensure(lines_.size(), cur_n_);
    if (ctx.greedy_bucket_next.size() != lines_.size()) ctx.greedy_bucket_next.resize(lines_.size());
    
    ctx.greedy_lines.clear();                    ctx.greedy_lines.reserve(lines_.size());
    ctx.candidate_choice.clear();                ctx.candidate_choice.reserve(lines_.size());
    ctx.disabled_lines.clear();                  ctx.disabled_lines.reserve(lines_.size());
    ctx.decremented_lines.clear();               ctx.decremented_lines.reserve(total_incidence_);
    ctx.newly_covered_points.clear();            ctx.newly_covered_points.reserve(as_size(cur_n_));
    ctx.removed_from_productive.clear();         ctx.removed_from_productive.reserve(lines_.size());
    ctx.removed_from_blocked_productive.clear(); ctx.removed_from_blocked_productive.reserve(lines_.size());
    ctx.state_blocked_ids.clear();               ctx.state_blocked_ids.reserve(64);
    ctx.saved_dual_points.clear();               ctx.saved_dual_points.reserve(as_size(cur_n_));
    ctx.saved_dual_values.clear();               ctx.saved_dual_values.reserve(as_size(cur_n_));
    
    ctx.forced_touched_lines.clear();            ctx.forced_touched_lines.reserve(lines_.size());
    ctx.covered_point_degree.clear();            ctx.covered_point_degree.reserve(as_size(cur_n_));
    ctx.covered_point_sole.clear();              ctx.covered_point_sole.reserve(as_size(cur_n_));
    ctx.include_ub_touched_lines.clear();        ctx.include_ub_touched_lines.reserve(256);
    ctx.undo_token_stack.clear();
    ctx.next_undo_token = 0;
    ctx.point_metadata_points.clear();           ctx.point_metadata_points.reserve(total_incidence_);
    ctx.point_metadata_degree.clear();           ctx.point_metadata_degree.reserve(total_incidence_);
    ctx.point_metadata_sole.clear();             ctx.point_metadata_sole.reserve(total_incidence_);
    ctx.point_metadata_prev_token.clear();       ctx.point_metadata_prev_token.reserve(total_incidence_);
    if (as_int(ctx.point_metadata_saved_token.size()) != cur_n_) {
        ctx.point_metadata_saved_token.resize(as_size(cur_n_));
    }
    std::fill(ctx.point_metadata_saved_token.begin(), ctx.point_metadata_saved_token.end(), -1);
    
    if (ctx.include_ub_line_reduction.size() != lines_.size()) ctx.include_ub_line_reduction.resize(lines_.size());
    if (initial_dual_seed_ready_) {
        const std::size_t seed_count = std::min<std::size_t>(initial_dual_seed_.size(), as_size(cur_n_));
        for (std::size_t i = 0; i < seed_count; ++i) {
            ctx.lagrangian.best_y[i] = initial_dual_seed_[i];
        }
        for (std::size_t i = seed_count; i < as_size(cur_n_); ++i) {
            ctx.lagrangian.best_y[i] = DualValue{0.0};
        }
        ctx.lagrangian.dual_ready = true;
    } else {
        ctx.lagrangian.dual_ready = false;
    }
}

void ExactGainSolver::rebuild_productive_point_metadata(Task& task) const {
    std::ranges::fill(task.productive_degree, 0);
    std::ranges::fill(task.sole_productive_line, kNoProductiveLine);
    
    for (const int id : task.productive_ids) {
        for (const CompactPointIndex point_cp : line_points(id)) {
            const int p = static_cast<int>(point_cp);
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;
            
            int& degree = task.productive_degree[as_size(p)];
            task.sole_productive_line[as_size(p)] = (degree == 0) ? id : kManyProductiveLines;
            ++degree;
        }
    }
}

void ExactGainSolver::save_point_metadata_for_undo(const Task& task, WorkerContext& ctx, int p) const {
    if (ctx.undo_token_stack.empty()) return;

    auto& saved_token = ctx.point_metadata_saved_token[as_size(p)];
    const int active_token = ctx.undo_token_stack.back();
    if (saved_token == active_token) return;

    ctx.point_metadata_points.push_back(p);
    ctx.point_metadata_degree.push_back(task.productive_degree[as_size(p)]);
    ctx.point_metadata_sole.push_back(task.sole_productive_line[as_size(p)]);
    ctx.point_metadata_prev_token.push_back(saved_token);
    saved_token = active_token;
}

#if PC_VALIDATE_FORCE_METADATA
bool ExactGainSolver::productive_point_metadata_matches_rebuild(const Task& task, int* mismatch_point) const {
    std::vector<int> expected_degree(as_size(cur_n_), 0);
    std::vector<int> expected_sole(as_size(cur_n_), kNoProductiveLine);

    for (const int id : task.productive_ids) {
        for (const CompactPointIndex point_cp : line_points(id)) {
            const int p = static_cast<int>(point_cp);
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;

            int& degree = expected_degree[as_size(p)];
            expected_sole[as_size(p)] = (degree == 0) ? id : kManyProductiveLines;
            ++degree;
        }
    }

    bool ok = true;
    task.active.for_each_set_bit(words_, [&](int p) {
        if (!ok) return;
        if (task.productive_degree[as_size(p)] != expected_degree[as_size(p)] ||
            task.sole_productive_line[as_size(p)] != expected_sole[as_size(p)]) {
            ok = false;
            if (mismatch_point) *mismatch_point = p;
        }
    });
    return ok;
}

void ExactGainSolver::probe_validate_metadata(const Task& task, const char* stage) const {
    int mismatch_point = -1;
    if (productive_point_metadata_matches_rebuild(task, &mismatch_point)) return;

    int expected_degree = 0;
    int expected_sole = kNoProductiveLine;
    for (const int id : task.productive_ids) {
        bool contains = false;
        for (const CompactPointIndex cp : line_points(id)) {
            const int p = static_cast<int>(cp);
            if (p >= cur_n_) break;
            if (p == mismatch_point) {
                contains = true;
                break;
            }
        }
        if (!contains) continue;
        expected_sole = (expected_degree == 0) ? id : kManyProductiveLines;
        ++expected_degree;
    }

    std::println(stderr,
                 "[probe] {} mismatch at N={} point={} active={} productive={} actual_degree={} expected_degree={} actual_sole={} expected_sole={} choice_size={}",
                 stage,
                 cur_n_,
                 mismatch_point,
                 task.active_count,
                 task.productive_ids.size(),
                 task.productive_degree[as_size(mismatch_point)],
                 expected_degree,
                 task.sole_productive_line[as_size(mismatch_point)],
                 expected_sole,
                 task.current_choice.size());
    std::abort();
}
#endif

int ExactGainSolver::find_any_productive_line_for_point(const Task& task, int p) const {
    // `ordered_point_lines_[p]` is pre-filtered by `activate_at <= cur_n_`.
    // Blocked lines or those with coverage < 3 are skipped via `productive_pos`.
    for (const int id : ordered_point_lines_[as_size(p)]) {
        if (task.productive_pos[as_size(id)] >= 0) return id;
    }
    return kNoProductiveLine;
}

void ExactGainSolver::on_productive_line_removed(Task& task, WorkerContext* ctx, int id) const {
    for (const CompactPointIndex point_cp : line_points(id)) {
        const int p = static_cast<int>(point_cp);
        if (p >= cur_n_) break;
        if (!task.active.test(p)) continue;

        int& degree = task.productive_degree[as_size(p)];
        int& sole   = task.sole_productive_line[as_size(p)];
        if (degree <= 0) continue;

        if (ctx) save_point_metadata_for_undo(task, *ctx, p);
        if (degree == 1) { 
            degree = 0; sole = kNoProductiveLine; 
        } else if (degree == 2) { 
            degree = 1; sole = find_any_productive_line_for_point(task, p); 
        } else { 
            --degree; sole = kManyProductiveLines; 
        }
    }
}

void ExactGainSolver::on_productive_line_reinserted(Task& task, WorkerContext* ctx, int id) const {
    for (const CompactPointIndex point_cp : line_points(id)) {
        const int p = static_cast<int>(point_cp);
        if (p >= cur_n_) break;
        if (!task.active.test(p)) continue;
        
        int& degree = task.productive_degree[as_size(p)];
        if (ctx) save_point_metadata_for_undo(task, *ctx, p);
        if (degree == 0) { 
            degree = 1; 
            task.sole_productive_line[as_size(p)] = id; 
        } else { 
            ++degree; 
            task.sole_productive_line[as_size(p)] = kManyProductiveLines; 
        }
    }
}

void ExactGainSolver::remove_productive_impl(Task& task, WorkerContext* ctx, int id, std::vector<int>* removed_log) const {
    int& pos = task.productive_pos[as_size(id)];
    if (pos < 0) return;
    
    const int last_id = task.productive_ids.back();
    task.productive_ids[as_size(pos)] = last_id;
    task.productive_pos[as_size(last_id)] = pos;
    task.productive_ids.pop_back();
    pos = -1;
    
    if (removed_log) removed_log->push_back(id);
    on_productive_line_removed(task, ctx, id);
}

void ExactGainSolver::remove_productive(Task& task, WorkerContext& ctx, int id) const { 
    remove_productive_impl(task, &ctx, id, &ctx.removed_from_productive); 
}

void ExactGainSolver::reinsert_productive(Task& task, WorkerContext* ctx, int id) const {
    if (task.productive_pos[as_size(id)] >= 0) return;
    
    task.productive_pos[as_size(id)] = as_int(task.productive_ids.size());
    task.productive_ids.push_back(id);
    on_productive_line_reinserted(task, ctx, id);
}

void ExactGainSolver::add_blocked_productive_impl(Task& task, int id) const {
    int& pos = task.blocked_productive_pos[as_size(id)];
    if (pos >= 0) return;
    
    auto& vec = task.blocked_productive_ids;
    auto it = std::lower_bound(vec.begin(), vec.end(), id);
    pos = static_cast<int>(it - vec.begin());
    vec.insert(it, id);
    
    for (std::size_t i = pos + 1; i < vec.size(); ++i) {
        task.blocked_productive_pos[as_size(vec[i])] = static_cast<int>(i);
    }
}

void ExactGainSolver::remove_blocked_productive_impl(Task& task, int id, std::vector<int>* removed_log) const {
    int& pos = task.blocked_productive_pos[as_size(id)];
    if (pos < 0) return;
    
    auto& vec = task.blocked_productive_ids;
    vec.erase(vec.begin() + pos);
    
    for (std::size_t i = pos; i < vec.size(); ++i) {
        task.blocked_productive_pos[as_size(vec[i])] = static_cast<int>(i);
    }
    pos = -1;
    if (removed_log) removed_log->push_back(id);
}

void ExactGainSolver::remove_blocked_productive(Task& task, WorkerContext& ctx, int id) const {
    remove_blocked_productive_impl(task, id, &ctx.removed_from_blocked_productive);
}

void ExactGainSolver::disable_line_without_undo(Task& task, int id) const {
    if (!task.available[as_size(id)]) return;
    task.available[as_size(id)] = 0;
    
    if (task.line_cover[as_size(id)] >= 3) {
        add_blocked_productive_impl(task, id);
    }
    remove_productive_impl(task, nullptr, id, nullptr);
}

void ExactGainSolver::disable_line(Task& task, WorkerContext& ctx, int id) const {
    if (!task.available[as_size(id)]) return;
    task.available[as_size(id)] = 0;
    
    if (task.line_cover[as_size(id)] >= 3) {
        add_blocked_productive_impl(task, id);
    }
    ctx.disabled_lines.push_back(id);
    remove_productive(task, ctx, id);
}

void ExactGainSolver::cover_point(Task& task, WorkerContext& ctx, int p) const {
    if (!task.active.test(p)) return;
    
    save_point_metadata_for_undo(task, ctx, p);
    task.productive_degree[as_size(p)]    = 0;
    task.sole_productive_line[as_size(p)] = kNoProductiveLine;
    task.active.reset(p);
    --task.active_count;
    ctx.newly_covered_points.push_back(p);
    
    for (const int id : ordered_point_lines_[as_size(p)]) {
        const bool available = task.available[as_size(id)] != 0;
        const bool blocked_productive = task.blocked_productive_pos[as_size(id)] >= 0;
        if (!available && !blocked_productive) continue;
        
        const uint16_t old_cover = task.line_cover[as_size(id)]--;
        ctx.decremented_lines.push_back(id);
        if (old_cover != 3) continue;
        
        if (available) remove_productive(task, ctx, id); 
        else remove_blocked_productive(task, ctx, id);
    }
}

void ExactGainSolver::select_line(Task& task, WorkerContext& ctx, int id) const {
    disable_line(task, ctx, id);
    for (const CompactPointIndex point_cp : line_points(id)) {
        const int p = static_cast<int>(point_cp);
        if (p >= cur_n_) break;
        cover_point(task, ctx, p);
    }
}

void ExactGainSolver::select_line_without_undo(Task& task, int id) const {
    if (!task.available[as_size(id)]) return;
    disable_line_without_undo(task, id);

    for (const CompactPointIndex point_cp : line_points(id)) {
        const int p = static_cast<int>(point_cp);
        if (p >= cur_n_) break;
        if (!task.active.test(p)) continue;

        task.active.reset(p);
        --task.active_count;
        
        for (const int lid : ordered_point_lines_[as_size(p)]) {
            const bool available = task.available[as_size(lid)] != 0;
            const bool blocked_productive = task.blocked_productive_pos[as_size(lid)] >= 0;
            if (!available && !blocked_productive) continue;
            if (task.line_cover[as_size(lid)]-- != 3) continue;
            
            if (available) remove_productive_impl(task, nullptr, lid, nullptr); 
            else remove_blocked_productive_impl(task, lid, nullptr);
        }
    }
}

void ExactGainSolver::undo_to(Task& task, WorkerContext& ctx, UndoMark mark) const {
    for (std::size_t i = mark.decremented_lines; i < ctx.decremented_lines.size(); ++i) {
        ++task.line_cover[as_size(ctx.decremented_lines[i])];
    }
    
    for (std::size_t i = mark.newly_covered_points; i < ctx.newly_covered_points.size(); ++i) {
        task.active.set(ctx.newly_covered_points[i]);
        ++task.active_count;
    }
    
    for (std::size_t i = mark.removed_from_blocked_productive; i < ctx.removed_from_blocked_productive.size(); ++i) {
        add_blocked_productive_impl(task, ctx.removed_from_blocked_productive[i]);
    }
        
    for (std::size_t i = mark.disabled_lines; i < ctx.disabled_lines.size(); ++i) {
        remove_blocked_productive_impl(task, ctx.disabled_lines[i], nullptr);
        task.available[as_size(ctx.disabled_lines[i])] = 1;
    }
    
    for (std::size_t i = mark.removed_from_productive; i < ctx.removed_from_productive.size(); ++i) {
        reinsert_productive(task, &ctx, ctx.removed_from_productive[i]);
    }

    for (std::size_t i = ctx.point_metadata_points.size(); i-- > mark.point_metadata_points;) {
        const int p = ctx.point_metadata_points[i];
        task.productive_degree[as_size(p)] = ctx.point_metadata_degree[i];
        task.sole_productive_line[as_size(p)] = ctx.point_metadata_sole[i];
        ctx.point_metadata_saved_token[as_size(p)] = ctx.point_metadata_prev_token[i];
    }

#if PC_VALIDATE_FORCE_METADATA
    probe_validate_metadata(task, "undo_to");
#endif
    ctx.removed_from_productive.resize(mark.removed_from_productive);
    ctx.removed_from_blocked_productive.resize(mark.removed_from_blocked_productive);
    ctx.decremented_lines.resize(mark.decremented_lines);
    ctx.newly_covered_points.resize(mark.newly_covered_points);
    ctx.covered_point_degree.resize(mark.covered_point_metadata);
    ctx.covered_point_sole.resize(mark.covered_point_metadata);
    ctx.point_metadata_points.resize(mark.point_metadata_points);
    ctx.point_metadata_degree.resize(mark.point_metadata_points);
    ctx.point_metadata_sole.resize(mark.point_metadata_points);
    ctx.point_metadata_prev_token.resize(mark.point_metadata_points);
    ctx.undo_token_stack.resize(mark.undo_token_depth_before);
    ctx.disabled_lines.resize(mark.disabled_lines);
}

ExactGainSolver::StateKey ExactGainSolver::make_state_key(const Task& task, WorkerContext& ctx) const {
    ctx.state_blocked_ids = task.blocked_productive_ids;
    return {task.active, ctx.state_blocked_ids};
}

bool ExactGainSolver::remember_frontier_state(const Task& task, WorkerContext& ctx, StateBestMap& best_gain) const {
    StateKey key = make_state_key(task, ctx);
    if (auto it = best_gain.find(key); it != best_gain.end()) {
        if (it->second >= task.current_gain) return true;
        it->second = task.current_gain;
        return false;
    }
    best_gain.emplace(std::move(key), task.current_gain);
    return false;
}

bool ExactGainSolver::frontier_state_is_stale(const Task& task, WorkerContext& ctx, const StateBestMap& best_gain) const {
    StateKey key = make_state_key(task, ctx);
    auto it = best_gain.find(key);
    return it != best_gain.end() && it->second > task.current_gain;
}

// ---------------------------------------------------------------------------
// Cover Inequality Lower Bound
// Finds the maximum set S ⊆ [n] such that no heavy line covers > 2 points of S.
// The value ceil(|S| / 2) provides a valid lower bound since any valid line
// (heavy, 2-point, or 1-point) can cover at most 2 points of S.
// This tightens the root LP relaxation more effectively than the Lagrangian alone.
//
// Algorithm: Greedily remove the point present in the most violated lines 
// (lines covering >= 3 points of S) until S is valid.
// ---------------------------------------------------------------------------
[[nodiscard]] int cover_inequality_lb(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
{
    if (n == 0) return 0;

    std::vector<unsigned char> in_support(as_size(n), 1);
    if (active_ids.empty()) return (n + 1) / 2;

    const int m = as_int(active_ids.size());
    std::vector<std::vector<int>> local_inc(as_size(n));
    std::vector<int> line_size(as_size(m), 0);
    
    for (int i = 0; i < m; ++i) {
        for (const int p : lines[as_size(active_ids[as_size(i)])].points) {
            if (p >= n) break;
            local_inc[as_size(p)].push_back(i);
            ++line_size[as_size(i)];
        }
    }

    int S_size = n;
    std::vector<int> violation_count(as_size(n), 0);
    for (int i = 0; i < m; ++i) {
        if (line_size[as_size(i)] >= 3) {
            for (const int p : lines[as_size(active_ids[as_size(i)])].points) {
                if (p >= n) break;
                ++violation_count[as_size(p)];
            }
        }
    }

    while (true) {
        int worst_p = -1, worst_v = 0;
        for (int p = 0; p < n; ++p) {
            if (in_support[as_size(p)] && violation_count[as_size(p)] > worst_v) {
                worst_v = violation_count[as_size(p)];
                worst_p = p;
            }
        }
        if (worst_v == 0) break;

        in_support[as_size(worst_p)] = 0;
        --S_size;

        for (const int i : local_inc[as_size(worst_p)]) {
            const bool was_violated = (line_size[as_size(i)] >= 3);
            if (--line_size[as_size(i)] < 3 && was_violated) {
                for (const int q : lines[as_size(active_ids[as_size(i)])].points) {
                    if (q >= n) break;
                    if (in_support[as_size(q)]) --violation_count[as_size(q)];
                }
            }
        }
    }

    return (S_size + 1) / 2;
}

[[nodiscard]] ExactSolveResult solve_exact(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines,
    const std::vector<std::vector<int>>& incidence, int root_lb_floor,
    const std::vector<int>& greedy_seed_line_ids, std::vector<int>& warm_start_solution,
    std::vector<double>& warm_dual_seed)
{
    ExactSolveResult out;
    const int cover_lb         = cover_inequality_lb(n, active_ids, lines);
    const int effective_lb_floor = std::max(root_lb_floor, cover_lb);
    
    ExactGainSolver solver(lines, incidence);
    GainSolveResult result     = solver.solve_for_n(n, active_ids, warm_start_solution, warm_dual_seed, greedy_seed_line_ids, effective_lb_floor);
    
    out.answer          = (n - result.best_gain + 1) / 2;
    out.lb              = result.root_cost_lb;
    out.lb_raw          = result.raw_root_cost_lb;
    out.lb_cov          = cover_lb;
    out.prod            = result.root_productive_lines;
    out.pinc            = result.root_productive_incidence;
    out.forced          = result.forced;
    out.forced_root     = result.forced_root;
    out.lag_iters       = result.lag_iters;
    out.polish_sweeps   = result.polish_sweeps;
    out.depth           = result.depth_max;
    out.lag_prune       = result.lag_prune;
    out.strong_branch   = result.strong_branch;
    out.nodes           = result.nodes;
    out.frontier_size   = result.frontier_size;
    out.root_ub_frac    = result.root_ub_frac;
    out.chosen_line_ids = std::move(result.chosen_line_ids);
    out.timed_out       = result.timed_out;
    
    warm_start_solution = out.chosen_line_ids;
    if (!result.dual_seed.empty()) {
        warm_dual_seed = std::move(result.dual_seed);
    }
    return out;
}

} // namespace solver

// ===========================================================================
// User Interface
// ===========================================================================
namespace ui {
// ---------------------------------------------------------------------------
// Background thread for real-time stderr progress reporting.
// Activates only for instances taking > 0.2 seconds to keep fast runs clean.
// Uses relaxed atomics (`begin_n`, `clear`) to communicate with the main thread 
// without penalizing the critical timing path. Outputs final results to stdout.
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
        using namespace std::chrono_literals;
        while (active_.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(100ms);
            if (!active_.load(std::memory_order_relaxed)) break;

            const int n = current_n_.load(std::memory_order_relaxed);
            if (n == 0) continue;

            const std::int64_t start_ns = n_start_ns_.load(std::memory_order_relaxed);
            if (start_ns == 0) continue;

            const std::int64_t now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - program_start_).count();
            const double n_elapsed    = static_cast<double>(now_ns - start_ns) * 1e-9;
            
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

// ── Reporter helpers ─────────────────────────────────────────────────────────
[[nodiscard]] int display_cols(const std::string& s) {
    int cols = 0;
    for (unsigned char c : s)
        if ((c & 0xC0) != 0x80) ++cols;
    return cols;
}

[[nodiscard]] std::string ordinal_suffix(int n) {
    const int mod100 = n % 100;
    if (mod100 >= 11 && mod100 <= 13) return "th";
    const int mod10 = n % 10;
    if (mod10 == 1) return "st";
    if (mod10 == 2) return "nd";
    if (mod10 == 3) return "rd";
    return "th";
}

[[nodiscard]] std::string format_commas(long long v) {
    if (v == 0) return "0";
    std::string s = std::to_string(v);
    int insert_pos = static_cast<int>(s.length()) - 3;
    while (insert_pos > (s[0] == '-' ? 1 : 0)) {
        s.insert(insert_pos, ",");
        insert_pos -= 3;
    }
    return s;
}

[[nodiscard]] std::string build_cover_string(
    int n,
    const std::vector<int>&             primes,
    const std::vector<core::HeavyLine>& lines,
    const std::vector<int>&             chosen_heavy_ids) {
    
    std::vector<bool> covered(util::as_size(n), false);
    std::string out = "[";
    bool first_line = true;
    
    for (int id : chosen_heavy_ids) {
        std::vector<int> current_line_points;
        for (int p : lines[util::as_size(id)].points) {
            if (p < n) {
                covered[util::as_size(p)] = true;
                current_line_points.push_back(p);
            }
        }
        std::ranges::sort(current_line_points);
        
        if (!first_line) out += ", ";
        first_line = false;
        
        out += "[";
        bool first_pt = true;
        for (int p : current_line_points) {
            if (!first_pt) out += ", ";
            out += "(" + std::to_string(p + 1) + ", " + std::to_string(primes[util::as_size(p)]) + ")";
            first_pt = false;
        }
        out += "]";
    }
    
    int pending_p = -1;
    for (int p = 0; p < n; ++p) {
        if (!covered[util::as_size(p)]) {
            if (pending_p == -1) {
                pending_p = p;
            } else {
                if (!first_line) out += ", ";
                first_line = false;
                out += "[(" + std::to_string(pending_p + 1) + ", " + std::to_string(primes[util::as_size(pending_p)]) + "), " +
                       "(" + std::to_string(p + 1) + ", " + std::to_string(primes[util::as_size(p)]) + ")]";
                pending_p = -1;
            }
        }
    }
    
    if (pending_p != -1) {
        if (!first_line) out += ", ";
        out += "[(" + std::to_string(pending_p + 1) + ", " + std::to_string(primes[util::as_size(pending_p)]) + ")]";
    }
    
    out += "]";
    return out;
}

} // namespace

// Entry point: parses args, generates geometry, and sweeps over N.
int main(int argc, char** argv) {
    using namespace config;
    using namespace util;
    using namespace core;
    using namespace solver;
    using namespace ui;

    const auto total_start = std::chrono::steady_clock::now();
    int requested_n = kExecutionLimit;
    if (argc >= 2) {
        if (const int parsed = std::atoi(argv[1]); parsed > 0) {
            requested_n = parsed;
        }
    }
    requested_n = std::min({requested_n, kExecutionLimit, kBitCapacity});

    const std::vector<int>       primes = generate_primes(requested_n);
    const std::vector<HeavyLine> lines  = enumerate_heavy_lines(primes);

    std::vector<std::vector<int>> incidence(as_size(requested_n));
    std::vector<std::vector<int>> activation(as_size(requested_n + 1));
    
    for (int id = 0; id < as_int(lines.size()); ++id) {
        for (const int p : lines[as_size(id)].points) {
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

    // ── Per-N record for supplementary summary ───────────────────────────────────
    struct NRecord {
        int       n;
        int       prime;
        int       answer;
        long long nodes;
        char      mode;
        double    elapsed;
        bool      is_awkward;
    };
    std::vector<NRecord> records;
    records.reserve(as_size(requested_n));

    int       prev                  = (kStartN > 1) ? kStartN - 1 : 0;
    long long total_nodes_processed = 0;
    bool      witness_ready         = false;
    
    std::vector<int> warm_start_solution;
    std::vector<double> warm_dual_seed;
    WitnessCover     witness_cover;
    
    for (int n_init = 1; n_init < kStartN; ++n_init) {
        for (const int id : activation[as_size(n_init)]) {
            active_ids.push_back(id);
        }
    }

    for (int n = std::max(kStartN, 1); n <= requested_n; ++n) {
        for (const int id : activation[as_size(n)]) {
            active_ids.push_back(id);
        }

        display.begin_n(n);
        const auto t0 = std::chrono::steady_clock::now();

        const GreedyCoverSolution greedy = greedy_upper_bound_solution(n, active_ids, lines);
        const int ub_raw                 = greedy.cost;
        int incumbent                    = std::min({n, prev + 1, greedy.cost});
        const int lb0_floor              = warm_start_solution.empty() ? 0 : prev;
        const std::size_t active_count   = active_ids.size();
        
        const bool witness_hit = witness_ready && witness_covers_point(witness_cover, n, primes[as_size(n - 1)]);
        if (witness_hit) incumbent = std::min(incumbent, prev);
        const int ub0_val = incumbent;

        // Capture warm-start line count before the witness check or solve_exact call.
        const int warm_size_val = as_int(warm_start_solution.size());
        ExactSolveResult exact_result;
        if (witness_hit) {
            exact_result.answer = prev;
            exact_result.lb     = prev;
            exact_result.lb_raw = prev;
        } else {
            exact_result = solve_exact(n, active_ids, lines, incidence, lb0_floor, greedy.line_ids, warm_start_solution, warm_dual_seed);
            if (!exact_result.timed_out) {
                witness_cover = build_witness_cover(n, requested_n, primes, lines, warm_start_solution);
                witness_ready = true;
            }
        }
        // Report the pre-solve warm-start size for every mode, including witness and root-only closures.
        exact_result.warm_size = warm_size_val;
        // Report how many heavy lines the greedy upper bound selected before any witness or exact solve decision.
        exact_result.greedy_heavy = as_int(greedy.line_ids.size());
        
        const int gap          = ub0_val - exact_result.lb;
        const int gap_raw      = ub_raw - exact_result.lb_raw;
        const char mode        = witness_hit ? 'W' : (exact_result.nodes == 0 ? 'R' : 'D');
        const double elapsed   = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        display.clear();
        total_nodes_processed += exact_result.nodes;

        const int forced              = mode == 'D' ? exact_result.forced : 0;
        const long long lag_iters     = mode == 'D' ? exact_result.lag_iters : 0;
        const int depth               = mode == 'D' ? exact_result.depth : 0;
        const int lag_prune           = mode == 'D' ? exact_result.lag_prune : 0;
        const int strong_branch       = mode == 'D' ? exact_result.strong_branch : 0;
        // New fields: four are D-only (solver internals); warm_size and greedy_heavy are valid for all modes.
        const int frontier_out        = mode == 'D' ? exact_result.frontier_size : 0;
        const long long polish_out    = mode == 'D' ? exact_result.polish_sweeps : 0LL;
        const double root_frac_out    = mode == 'D' ? exact_result.root_ub_frac : 0.0;
        const int forced_root_out     = mode == 'D' ? exact_result.forced_root : 0;

        if (exact_result.timed_out) {
            std::println("[stopped] N={} time limit of {}s exceeded (cancelled mid-run after {})", 
                         n, kPerNTimeLimitSeconds, format_seconds(elapsed));
            break;
        }

        // ── Supplementary human-readable output (OEIS format) ────────────────────────
        bool is_awkward = (exact_result.answer > prev) && !records.empty();
        records.push_back({n, primes[as_size(n - 1)], exact_result.answer, exact_result.nodes, mode, elapsed, is_awkward});

        const double total_so_far = std::chrono::duration<double>(std::chrono::steady_clock::now() - program_start).count();
        std::string core_str = "── " + std::to_string(n) + ordinal_suffix(n) + 
                               " prime (" + std::to_string(primes[as_size(n - 1)]) + ") · " + 
                               std::to_string(exact_result.answer) + " lines · " + 
                               format_seconds(elapsed) + " · total: " +
                               format_seconds(total_so_far) + " ";
        const int pad_cols = std::max(0, 80 - display_cols(core_str));
        for (int i = 0; i < pad_cols; ++i) core_str += "─";
        std::println("{}", core_str);

        std::println("   Stats:  N={} prime={} lines={} time={} mode={} active={} nodes={} frontier={}"
                     " ub0={} lb_floor={} lb_cov={} lb={} gap={} ub_raw={} lb_raw={} gap_raw={}"
                     " root_ub_frac={:.3f} warm_size={} greedy_heavy={} prod={} pinc={}"
                     " forced_root={} forced={} lag_iters={} polish_sweeps={} lag_prune={}"
                     " strong_branch={} depth={}",
            n,
            primes[as_size(n - 1)],
            exact_result.answer,
            format_seconds(elapsed),
            mode,
            active_count,
            exact_result.nodes,
            frontier_out,
            ub0_val,
            lb0_floor,
            exact_result.lb_cov,
            exact_result.lb,
            gap,
            ub_raw,
            exact_result.lb_raw,
            gap_raw,
            root_frac_out,
            exact_result.warm_size,
            exact_result.greedy_heavy,
            exact_result.prod,
            exact_result.pinc,
            forced_root_out,
            forced,
            lag_iters,
            polish_out,
            lag_prune,
            strong_branch,
            depth);

        std::string cover_str = build_cover_string(n, primes, lines, 
            (mode == 'W') ? warm_start_solution : exact_result.chosen_line_ids);
        std::println("   Cover: {}", cover_str);
        std::println("");
        
        prev = exact_result.answer;

        if (kPerNTimeLimitSeconds > 0.0 && elapsed > kPerNTimeLimitSeconds) {
            std::println("[stopped] per-N time limit of {}s exceeded", kPerNTimeLimitSeconds);
            break;
        }
    }
    const double total_elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - total_start).count();
    std::println("total time={} total nodes={}", format_seconds(total_elapsed), total_nodes_processed);
    
    // ── Supplementary global summary ─────────────────────────────────────────────
    if (!records.empty()) {
        int first_n = records.front().n;
        int last_n = records.back().n;
        int count = static_cast<int>(records.size());
        
        int min_lines = records.front().answer;
        int min_n = records.front().n;
        int max_lines = records.front().answer;
        int max_n = records.front().n;
        
        int awk_count = 0;
        int max_n_chars = 0, max_p_chars = 0, max_l_chars = 0;
        for (const auto& r : records) {
            if (r.answer < min_lines) { min_lines = r.answer; min_n = r.n; }
            if (r.answer > max_lines) { max_lines = r.answer; max_n = r.n; }
            if (r.is_awkward) {
                ++awk_count;
                max_n_chars = std::max(max_n_chars, static_cast<int>(std::to_string(r.n).length()));
                max_p_chars = std::max(max_p_chars, static_cast<int>(std::to_string(r.prime).length()));
                max_l_chars = std::max(max_l_chars, static_cast<int>(std::to_string(r.answer).length()));
            }
        }
        
        std::println("================================================================================");
        std::println("PRIME COVER SEQUENCE  N={}..{}", first_n, last_n);
        std::println("Total time  : {}", format_seconds(total_elapsed));
        std::println("Total nodes : {}", format_commas(total_nodes_processed));
        std::println("N range     : {}..{}  ({} values)", first_n, last_n, count);
        std::println("Line counts : min={} (first at N={})  max={} (first at N={})", min_lines, min_n, max_lines, max_n);
        std::println("Awkward primes ({} total — N where line count increased from previous):", awk_count);
        
        if (awk_count == 0) {
            std::println("  none found in this run");
        } else {
            for (const auto& r : records) {
                if (r.is_awkward) {
                    std::string n_str = std::to_string(r.n);
                    std::string p_str = std::to_string(r.prime);
                    std::string l_str = std::to_string(r.answer);
                    std::println("  N={:<{}}  prime={:<{}}  lines={:<{}}  (+1)", 
                        n_str, max_n_chars, p_str, max_p_chars, l_str, max_l_chars);
                }
            }
        }
        
        if (last_n < requested_n) {
            std::println("================ (run stopped early — last completed N={}) ================", last_n);
        } else {
            std::println("================================================================================");
        }
    }
    
    return 0;
}