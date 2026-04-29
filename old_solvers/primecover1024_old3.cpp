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
 *   - MRV branching on the uncovered point with the fewest productive lines.
 *   - Forced propagation for points not lying on any productive heavy line.
 *   - Incremental state tracking: productive_count / size_freq updated via
 *     delta application instead of full recomputation at each DFS node.
 *   - Optimistic lower bound + branch-and-bound pruning.
 *   - Memoization on (active mask, dangling parity) in the monolithic solver.
 *   - Exact component fallback when the active heavy-line graph truly
 *     decomposes; dense rigid instances stay on the faster monolithic path.
 */

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr int kBitCapacity = 1024;
constexpr int kBitWords    = kBitCapacity / 64;
constexpr int kExecutionLimit = 159;
constexpr int kDefaultRunN    = 159;

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
// All fields are updated via apply_delta rather than recomputed from scratch
// at every node, which is the key source of the ~5x speedup over the
// full-recompute approach.
struct IncrementalState {
    std::vector<int> line_cover;       // current active intersection size per line
    std::vector<int> productive_count; // # productive lines covering each point
    std::vector<int> best_cover;       // largest productive line covering each point
    std::vector<int> size_freq;        // size_freq[k] = # productive lines with cover == k
};

struct ActiveComponent {
    std::vector<int> line_ids;
    BitMask1024      point_mask;
    int              point_count = 0;
};

using Summary2 = std::array<int, 2>;   // {even-parity cost, odd-parity cost}

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

[[nodiscard]] int greedy_upper_bound(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
{
    const int words = words_for_n(n);
    BitMask1024 active = BitMask1024::prefix(n);
    int count = n, cost = 0;
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
        ++cost;
    }
    return cost + ceil_half(count);
}

// ---------------------------------------------------------------------------
// MonolithicExactSolver
// Incremental state: productive_count / size_freq are maintained via
// apply_delta rather than recomputed at every DFS node.
// ---------------------------------------------------------------------------
class MonolithicExactSolver {
public:
    MonolithicExactSolver(const std::vector<HeavyLine>&            lines,
                          const std::vector<std::vector<int>>&     incidence)
        : lines_(lines), incidence_(incidence) {}

    [[nodiscard]] int solve_for_n(int n, const std::vector<int>& active_ids, int incumbent) {
        cur_n_  = n;
        words_  = words_for_n(n);
        best_   = incumbent;
        StateKeyWithParityHash::used_words = words_;
        memo_.clear();

        // Build the initial IncrementalState for the full active set.
        IncrementalState state;
        state.line_cover.assign(lines_.size(), 0);
        state.productive_count.assign(static_cast<std::size_t>(n), 0);
        state.best_cover.assign(static_cast<std::size_t>(n), 0);
        state.size_freq.assign(static_cast<std::size_t>(n + 1), 0);

        const BitMask1024 active = BitMask1024::prefix(n);
        for (int id : active_ids) {
            const int c = lines_[id].mask.intersection_count(active, words_);
            state.line_cover[static_cast<std::size_t>(id)] = c;
            if (c >= 3) {
                state.size_freq[static_cast<std::size_t>(c)]++;
                for (int p : lines_[id].points) if (p < n) {
                    state.productive_count[static_cast<std::size_t>(p)]++;
                    state.best_cover[static_cast<std::size_t>(p)] =
                        std::max(state.best_cover[static_cast<std::size_t>(p)], c);
                }
            }
        }

        dfs(active, n, false, 0, std::move(state));
        return best_;
    }

private:
    // Propagate the removal of all bits set in `removed` into state `s`.
    // For each removed point, decrement the cover of every line covering it.
    // When a line's cover drops below 3 it is no longer productive, so
    // decrement productive_count for its remaining points.
    void apply_delta(IncrementalState& s, const BitMask1024& removed) const {
        removed.for_each_set_bit(words_, [&](int p) {
            for (int id : incidence_[static_cast<std::size_t>(p)]) {
                if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
                int& lc = s.line_cover[static_cast<std::size_t>(id)];
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

    // Fractional lower bound: greedily pack the largest available lines and
    // cost remaining isolated points at ceil(residual/2).
    [[nodiscard]] int lower_bound(const IncrementalState& s, int ac, bool d) const {
        const int pending = ac + static_cast<int>(d);
        int best = (pending + 1) / 2, covered = 0, chosen = 0, max_c = 2;
        for (int i = ac; i >= 3; --i)
            if (s.size_freq[static_cast<std::size_t>(i)] > 0) { max_c = i; break; }
        for (int i = ac; i >= 3; --i) {
            for (int j = 0; j < s.size_freq[static_cast<std::size_t>(i)]; ++j) {
                ++chosen;
                covered = std::min(ac, covered + i);
                best    = std::min(best, chosen + ceil_half(std::max(0, pending - covered)));
                if (covered == ac) return std::max(best, (pending + max_c - 1) / max_c);
            }
        }
        return std::max(best, (pending + max_c - 1) / max_c);
    }

    void dfs(BitMask1024 active, int ac, bool d, int cost, IncrementalState state) {
        best_ = std::min(best_, cost + (ac + static_cast<int>(d) + 1) / 2);
        if (cost >= best_) return;

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
                best_ = std::min(best_, cost + (ac + static_cast<int>(d) + 1) / 2);
                return;
            }
            if (fc == 0) {
                if (cost + lower_bound(state, ac, d) >= best_) return;
                break;
            }
            apply_delta(state, forced);
            active = active.and_not(forced);
            ac -= fc;
            const int total = static_cast<int>(d) + fc;
            cost += total / 2;
            d     = (total & 1) != 0;
            if (cost >= best_) return;
        }
        if (ac == 0) { best_ = std::min(best_, cost + static_cast<int>(d)); return; }

        // Memoize on (active mask, dangling parity): prune revisits at equal
        // or higher cost.
        const StateKeyWithParity key{active.words, d};
        const auto [it, inserted] = memo_.emplace(key, cost);
        if (!inserted) {
            if (it->second <= cost) return;
            it->second = cost;
        }

        // MRV: branch on the point covered by the fewest productive lines.
        int bp = -1, mc = std::numeric_limits<int>::max();
        active.for_each_set_bit(words_, [&](int p) {
            const int cnt = state.productive_count[static_cast<std::size_t>(p)];
            if (cnt < mc) { mc = cnt; bp = p; }
        });

        // Enumerate lines covering bp, sort by descending current cover size.
        std::vector<int> candidates;
        for (int id : incidence_[static_cast<std::size_t>(bp)]) {
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
            if (state.line_cover[static_cast<std::size_t>(id)] >= 3)
                candidates.push_back(id);
        }
        std::sort(candidates.begin(), candidates.end(), [&](int a, int b) {
            return state.line_cover[static_cast<std::size_t>(a)] >
                   state.line_cover[static_cast<std::size_t>(b)];
        });

        for (int id : candidates) {
            const int             cover = state.line_cover[static_cast<std::size_t>(id)];
            const BitMask1024     rem   = active & lines_[static_cast<std::size_t>(id)].mask;
            IncrementalState      ns    = state;
            apply_delta(ns, rem);
            dfs(active.and_not(lines_[static_cast<std::size_t>(id)].mask),
                ac - cover, d, cost + 1, std::move(ns));
        }

        // Skip branch: defer bp as an isolated point.
        BitMask1024      skp{};
        skp.set(bp);
        IncrementalState skip_s = state;
        apply_delta(skip_s, skp);
        const int total = static_cast<int>(d) + 1;
        dfs(active.and_not(skp), ac - 1,
            (total & 1) != 0, cost + total / 2, std::move(skip_s));
    }

    const std::vector<HeavyLine>&        lines_;
    const std::vector<std::vector<int>>& incidence_;
    int  cur_n_ = 0, words_ = 0, best_ = 0;
    std::unordered_map<StateKeyWithParity, int, StateKeyWithParityHash> memo_;
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

[[nodiscard]] ComponentBuildResult build_active_components(
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
    const std::vector<std::vector<int>>& incidence, int incumbent)
{
    const ComponentBuildResult built = build_active_components(n, active_ids, lines);

    // Single component (or no lines): go straight to the monolithic solver.
    if (built.components.size() <= 1) {
        MonolithicExactSolver solver(lines, incidence);
        return solver.solve_for_n(n, active_ids, incumbent);
    }

    // Dense rigid webs: monolithic memoized DFS is faster than component DP.
    int largest = 0;
    for (const ActiveComponent& c : built.components)
        largest = std::max(largest, c.point_count);
    if (largest * 10 >= n * 9) {
        MonolithicExactSolver solver(lines, incidence);
        return solver.solve_for_n(n, active_ids, incumbent);
    }

    // Genuinely decomposed: solve each component independently and combine
    // via a parity-aware DP to track the global dangling count.
    const int inf = incumbent + 1;
    Summary2  dp{0, inf};

    for (const ActiveComponent& comp : built.components) {
        ComponentProfileSolver solver(lines, n, comp.line_ids, incumbent);
        const Summary2 summary = solver.solve(comp.point_mask);

        Summary2 next{inf, inf};
        next[0] = std::min(
            capped_add(dp[0], summary[0], incumbent),
            capped_add(capped_add(dp[1], summary[1], incumbent), 1, incumbent));
        next[1] = std::min(
            capped_add(dp[0], summary[1], incumbent),
            capped_add(dp[1], summary[0], incumbent));
        dp = next;
    }

    const Summary2 iso{
        ceil_half(built.isolated_points),
        built.isolated_points == 0 ? inf : built.isolated_points / 2
    };

    Summary2 final_dp{inf, inf};
    final_dp[0] = std::min(
        capped_add(dp[0], iso[0], incumbent),
        capped_add(capped_add(dp[1], iso[1], incumbent), 1, incumbent));
    final_dp[1] = std::min(
        capped_add(dp[0], iso[1], incumbent),
        capped_add(dp[1], iso[0], incumbent));

    return std::min(final_dp[0], capped_add(final_dp[1], 1, incumbent));
}

} // namespace

int main(int argc, char** argv) {
    int requested_n = kDefaultRunN;
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

    int prev = 0;
    for (int n = 1; n <= requested_n; ++n) {
        for (int id : activation[static_cast<std::size_t>(n)])
            active_ids.push_back(id);

        const auto t0 = std::chrono::steady_clock::now();

        int incumbent = std::min(n, prev + 1);
        incumbent     = std::min(incumbent, greedy_upper_bound(n, active_ids, lines));

        const int answer = solve_exact_with_components(
            n, active_ids, lines, incidence, incumbent);

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();

        std::cout << "N="     << n
                  << " prime=" << primes[static_cast<std::size_t>(n - 1)]
                  << " lines=" << answer
                  << " time="  << format_seconds(elapsed)
                  << '\n';
        prev = answer;
    }
    return 0;
}