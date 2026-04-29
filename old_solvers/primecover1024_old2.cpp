/*
 * primecover1024.cpp — C++26, Bit-Parallel Maximum Net Gain Solver
 * with Incremental Event-Driven Architecture and Connected-Component Decomposition
 *
 * HYBRID v4 — merges MRV (Fail-First) branching + LP-Dual Packing Bound from
 * the direct set-cover reformulation into the original incremental Net-Gain
 * architecture.  No essential features were removed.
 *
 *   ① MRV (Fail-First) Point-Centric Branching [NEW]:
 *       At each DFS node, scan all uncovered points and pick the one with the
 *       fewest productive covering lines (most constrained).  Replaces the
 *       original fixed-order linear scan; forces rapid tree collapse.
 *
 *   ② LP-Dual Packing Bound [NEW]:
 *       While computing MRV we accumulate dual_sum = Σ 1/max_sz(p) over
 *       uncovered points.  This gives a fractional lower bound on remaining
 *       lines, and hence an upper bound on remaining gain:
 *           remaining_gain ≤ uncov_pts − 2·⌈dual_sum⌉
 *       Prune when cur_gain + remaining_gain_ub ≤ best_gain.
 *
 *   ③ Size-sorted line branches tried BEFORE the skip branch [NEW]:
 *       Largest-overlap lines are explored first, tightening best_gain faster.
 *       (Original explored skip branch first.)
 *
 *   ④ N-1 Warm-Start forwarding [NEW]:
 *       cached_G_star_ from the previous N feeds the initial lower bound for
 *       the first component at N+1, restricting search depth from the start.
 *
 *   ⑤ Rollback-based DFS [from v3]:
 *       covered |= new_bits → recurse → covered.andnot_inplace(new_bits)
 *       covered.set(p)      → recurse → covered.clear_bit(p)
 *       Eliminates ALL Bitset1024 pass-by-value copies from the hot path.
 *
 *   ⑥ Greedy warm-start bug fix + size-sorted order [from v3]:
 *       comp pre-sorted by decreasing line size → tighter initial bound.
 *
 *   ⑦ Union-Find component decomposition [from v3]:
 *       Path-halving union-find replaces adj vector<vector<int>> + BFS.
 *
 *   ⑧ Live UI thread [from original]:
 *       Background thread prints current N and elapsed time at 100 ms cadence.
 *
 * Compile:
 *   g++ -std=c++26 -O3 -march=native -o primecover1024 primecover1024.cpp
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <span>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

#if __has_include(<print>)
#  include <print>
#else
#  include <format>
#  include <iostream>
namespace std {
    template<class... Args>
    void print(std::format_string<Args...> fmt, Args&&... args) {
        std::cout << std::format(fmt, std::forward<Args>(args)...);
    }
    template<class... Args>
    void println(std::format_string<Args...> fmt, Args&&... args) {
        std::cout << std::format(fmt, std::forward<Args>(args)...) << '\n';
    }
}
#endif

// ═══════════════════════════════════════════════════════════════════════════
// §0  Enhanced 1024-bit bitset with rollback primitives
// ═══════════════════════════════════════════════════════════════════════════

struct Bitset1024 {
    std::array<uint64_t, 16> words{};

    void set(int bit) noexcept {
        words[bit >> 6] |= (1ULL << (bit & 63));
    }

    // clear a single bit (used for skip-branch rollback)
    void clear_bit(int bit) noexcept {
        words[bit >> 6] &= ~(1ULL << (bit & 63));
    }

    bool test(int bit) const noexcept {
        return (words[bit >> 6] & (1ULL << (bit & 63))) != 0;
    }

    int count() const noexcept {
        int c = 0;
        for (uint64_t w : words) c += std::popcount(w);
        return c;
    }

    bool none() const noexcept {
        for (uint64_t w : words) if (w) return false;
        return true;
    }

    bool intersects(const Bitset1024& o) const noexcept {
        for (int i = 0; i < 16; ++i)
            if (words[i] & o.words[i]) return true;
        return false;
    }

    Bitset1024 operator|(const Bitset1024& o) const noexcept {
        Bitset1024 r;
        for (int i = 0; i < 16; ++i) r.words[i] = words[i] | o.words[i];
        return r;
    }

    Bitset1024 operator&(const Bitset1024& o) const noexcept {
        Bitset1024 r;
        for (int i = 0; i < 16; ++i) r.words[i] = words[i] & o.words[i];
        return r;
    }

    [[nodiscard]] Bitset1024 andnot(const Bitset1024& o) const noexcept {
        Bitset1024 r;
        for (int i = 0; i < 16; ++i) r.words[i] = words[i] & ~o.words[i];
        return r;
    }

    Bitset1024& operator|=(const Bitset1024& o) noexcept {
        for (int i = 0; i < 16; ++i) words[i] |= o.words[i];
        return *this;
    }

    // in-place andnot (used for line-branch rollback)
    void andnot_inplace(const Bitset1024& o) noexcept {
        for (int i = 0; i < 16; ++i) words[i] &= ~o.words[i];
    }

    Bitset1024 operator~() const noexcept {
        Bitset1024 r;
        for (int i = 0; i < 16; ++i) r.words[i] = ~words[i];
        return r;
    }
};

// ═══════════════════════════════════════════════════════════════════════════
// Utility
// ═══════════════════════════════════════════════════════════════════════════

[[nodiscard]]
static std::string fmt_dur(double seconds)
{
    if (seconds < 1e-3)  return std::format("{:.2f} µs",  seconds * 1e6);
    if (seconds < 1.0)   return std::format("{:.2f} ms",  seconds * 1e3);
    if (seconds < 60.0)  return std::format("{:.3f} s",   seconds);
    if (seconds < 3600.) return std::format("{:.2f} min", seconds / 60.);
    return                      std::format("{:.2f} hr",  seconds / 3600.);
}

// ═══════════════════════════════════════════════════════════════════════════
// §1  Segmented sieve
// ═══════════════════════════════════════════════════════════════════════════

[[nodiscard]]
static std::vector<int> generate_primes(int count)
{
    const double fn  = static_cast<double>(count);
    const int    lim = (count < 10)
        ? 50
        : static_cast<int>(fn * (std::log(fn) + std::log(std::log(fn)) + 2.5)) + 500;

    const int sq = static_cast<int>(std::sqrt(static_cast<double>(lim))) + 2;

    std::vector<bool> small(static_cast<std::size_t>(sq + 1), true);
    small[0] = small[1] = false;
    for (int p = 2; p * p <= sq; ++p)
        if (small[static_cast<std::size_t>(p)])
            for (int j = p * p; j <= sq; j += p)
                small[static_cast<std::size_t>(j)] = false;

    std::vector<int> base, primes;
    primes.reserve(static_cast<std::size_t>(count + 64));
    for (int p = 2; p <= sq; ++p)
        if (small[static_cast<std::size_t>(p)]) { base.push_back(p); primes.push_back(p); }

    constexpr int BLOCK = 1 << 15;
    for (int lo = sq + 1; static_cast<int>(primes.size()) < count; lo += BLOCK) {
        const int hi = std::min(lo + BLOCK - 1, lim);
        std::vector<bool> seg(static_cast<std::size_t>(hi - lo + 1), true);
        for (int p : base) {
            int start = static_cast<int>(
                (static_cast<long long>(lo) + p - 1) / p * p);
            if (start == p) start += p;
            for (int j = start; j <= hi; j += p)
                seg[static_cast<std::size_t>(j - lo)] = false;
        }
        for (int i = lo; i <= hi && static_cast<int>(primes.size()) < count; ++i)
            if (seg[static_cast<std::size_t>(i - lo)])
                primes.push_back(i);
    }
    primes.resize(static_cast<std::size_t>(count));
    return primes;
}

// ═══════════════════════════════════════════════════════════════════════════
// §2  Collinear line enumeration
// ═══════════════════════════════════════════════════════════════════════════

struct CollinearLine {
    std::vector<int> pts;         // point indices, sorted ascending
    int              activate_at; // pts[2]+1  (first n where line has ≥3 pts)
};

[[nodiscard]]
static std::vector<CollinearLine> enumerate_collinear_lines(
    const std::vector<int>& primes)
{
    const int N = static_cast<int>(primes.size());
    std::vector<Bitset1024> used(static_cast<std::size_t>(N));

    auto pack_slope = [](long long dy_r, long long dx_r) noexcept -> int {
        return static_cast<int>(dy_r * 1024 + dx_r);
    };

    std::vector<CollinearLine> all_lines;

    for (int i = 0; i < N; ++i) {
        std::unordered_map<int, std::vector<int>> slope_map;
        slope_map.reserve(static_cast<std::size_t>(N - i));

        for (int j = i + 1; j < N; ++j) {
            if (used[static_cast<std::size_t>(i)].test(j)) continue;
            const long long dy = static_cast<long long>(primes[j]) - primes[i];
            const long long dx = static_cast<long long>(j - i);
            const long long g  = std::gcd(dy, dx);
            slope_map[pack_slope(dy / g, dx / g)].push_back(j);
        }

        for (auto& [key, js] : slope_map) {
            if (static_cast<int>(js.size()) < 2) continue;

            std::vector<int> pts;
            pts.reserve(js.size() + 1);
            pts.push_back(i);
            pts.insert(pts.end(), js.begin(), js.end());

            for (int a = 0; a < static_cast<int>(pts.size()); ++a)
                for (int b = a + 1; b < static_cast<int>(pts.size()); ++b)
                    used[static_cast<std::size_t>(pts[a])].set(pts[b]);

            const int act = pts[2] + 1;
            all_lines.push_back({std::move(pts), act});
        }
    }
    return all_lines;
}

// ═══════════════════════════════════════════════════════════════════════════
// Line structure
// ═══════════════════════════════════════════════════════════════════════════

struct Line {
    Bitset1024 mask;
    int size;
    int max_pt;
    int id;
};

// ═══════════════════════════════════════════════════════════════════════════
// §3  Maximum Net-Gain component solver
//
//  HYBRID: replaces the fixed-order PointCentricDFS with MRV (Fail-First)
//  branching and an LP-Dual Packing Bound, while retaining:
//    • rollback-based DFS (no Bitset1024 pass-by-value)
//    • greedy warm-start with size-sorted component order
//    • point_potential LP-dual ceiling (original bound, kept as primary gate)
//    • packing bound as a secondary, tighter cut
// ═══════════════════════════════════════════════════════════════════════════

[[nodiscard]]
static int solve_max_gain_component(
    std::span<const std::vector<int>>  line_pts,
    std::span<const Bitset1024>        line_masks,
    std::vector<int>                   comp,
    int                                initial_lb = 0)
{
    const int c = static_cast<int>(comp.size());
    if (c == 0) return 0;
    if (c == 1) {
        return static_cast<int>(line_pts[comp[0]].size()) - 2;
    }

    // ── Sort comp by decreasing line size for better greedy ordering ────────
    std::ranges::sort(comp, [&](int a, int b) {
        return line_pts[a].size() > line_pts[b].size();
    });

    // ── Greedy warm-start (BUG FIX: use comp[i], not i; size-sorted above) ─
    int best_gain = initial_lb;
    {
        Bitset1024 gc{};
        int gw = 0;
        for (int ci : comp) {
            const int np = line_masks[ci].andnot(gc).count();
            if (np >= 3) {
                gc |= line_masks[ci];
                gw += np - 2;
            }
        }
        if (gw > best_gain) best_gain = gw;
    }

    // ── Build per-component data structures ────────────────────────────────
    std::vector<Line>              active_lines;
    active_lines.reserve(c);
    std::vector<std::vector<int>>  incidence_map(1024);
    std::vector<int>               point_potential(1024, 0);

    int local_id = 0;
    for (int ci : comp) {
        const auto& pts_list = line_pts[ci];
        const int sz = static_cast<int>(pts_list.size());
        Bitset1024 mask{};
        int max_ptt = -1;
        for (int pp : pts_list) {
            mask.set(pp);
            if (pp > max_ptt) max_ptt = pp;
        }
        // Max potential contribution per point: (sz-2)/sz * 60, scaled integer
        const int max_d = ((sz - 2) * 60) / sz;
        for (int pp : pts_list) {
            incidence_map[pp].push_back(local_id);
            if (max_d > point_potential[pp])
                point_potential[pp] = max_d;
        }
        active_lines.push_back({mask, sz, max_ptt, local_id});
        ++local_id;
    }

    // Sort incidence lists by decreasing line size (largest line tried first)
    for (int p = 0; p < 1024; ++p) {
        auto& inc = incidence_map[p];
        if (inc.size() > 1) {
            std::ranges::sort(inc, [&](int idA, int idB) {
                return active_lines[idA].size > active_lines[idB].size;
            });
        }
    }

    // Build the set of all active points (appear in ≥1 active line)
    std::vector<int> active_pts_vec;
    active_pts_vec.reserve(64);
    for (int p = 0; p < 1024; ++p)
        if (!incidence_map[p].empty()) active_pts_vec.push_back(p);

    if (active_pts_vec.empty()) return best_gain;

    int total_potential = 0;
    for (int p : active_pts_vec) total_potential += point_potential[p];

    // ── MRV Rollback-based DFS ─────────────────────────────────────────────
    //
    //  At each node, scan all uncovered points to find the "most constrained"
    //  one (fewest productive lines through it).  This is the MRV heuristic:
    //  branching on the most constrained variable first minimises tree size.
    //
    //  While scanning, accumulate dual_sum = Σ 1/max_sz(p), which gives a
    //  fractional lower bound on remaining lines and hence an upper bound on
    //  remaining gain:  remaining_gain ≤ uncov_pts − 2·⌈dual_sum⌉.
    //
    //  Two pruning gates, applied in order:
    //    Gate 1 (fast):   cur_gain + cur_pot/60 ≤ best_gain  [original bound]
    //    Gate 2 (tight):  cur_gain + uncov − 2·⌈dual_sum⌉ ≤ best_gain  [new]
    //
    //  Rollback is unchanged:
    //    line branch:  covered |= nb → recurse → covered.andnot_inplace(nb)
    //    skip branch:  covered.set(p) → recurse → covered.clear_bit(p)
    //
    struct PointCentricDFS {
        int&                                  best_gain;
        const std::vector<int>&               active_pts_vec;
        const std::vector<std::vector<int>>&  incidence_map;
        const std::vector<Line>&              active_lines;
        const std::vector<int>&               point_potential;
        Bitset1024                            covered{};   // ← mutable rollback state

        PointCentricDFS(int& bg,
                        const std::vector<int>& apv,
                        const std::vector<std::vector<int>>& im,
                        const std::vector<Line>& al,
                        const std::vector<int>& ppot)
            : best_gain(bg), active_pts_vec(apv)
            , incidence_map(im), active_lines(al), point_potential(ppot)
        {}

        void operator()(int cur_gain, int cur_pot) {
            // ── Gate 1: original LP-dual potential bound (fast check) ────
            if (cur_gain + (cur_pot / 60) <= best_gain) return;

            // ── MRV scan: find most constrained uncovered point ──────────
            //    Also accumulate dual_sum for Gate 2.
            int    best_p    = -1;
            int    min_cands = 1'000'000;
            int    best_p_pp = 0;
            double dual_sum  = 0.0;
            int    uncov_pts = 0;

            for (int p : active_pts_vec) {
                if (covered.test(p)) continue;
                ++uncov_pts;

                int cand_count = 0;
                int max_sz     = 2; // fallback: point will be covered as part of a pair

                for (int lid : incidence_map[p]) {
                    // Count uncovered points in this line's mask
                    const Bitset1024& lm = active_lines[lid].mask;
                    int sz = 0;
                    for (int wi = 0; wi < 16; ++wi)
                        sz += std::popcount(lm.words[wi] & ~covered.words[wi]);
                    if (sz >= 3) {
                        ++cand_count;
                        if (sz > max_sz) max_sz = sz;
                    }
                }
                // Fractional contribution of point p to lines lower bound
                dual_sum += 1.0 / static_cast<double>(max_sz);

                if (cand_count < min_cands) {
                    min_cands = cand_count;
                    best_p    = p;
                    best_p_pp = point_potential[p];
                }
            }

            // All points covered → update best
            if (best_p == -1) {
                if (cur_gain > best_gain) best_gain = cur_gain;
                return;
            }

            // ── Gate 2: LP-Dual Packing Bound (new, tighter) ────────────
            //   remaining lines   ≥ ⌈dual_sum⌉
            //   remaining gain    ≤ uncov_pts − 2·⌈dual_sum⌉
            {
                const int lines_lb   = static_cast<int>(std::ceil(dual_sum - 1e-9));
                const int gain_ub    = uncov_pts - 2 * lines_lb;
                if (cur_gain + gain_ub <= best_gain) return;
            }

            // ── Branch 1…k: try each productive line through best_p ──────
            //    Incidence list is pre-sorted by decreasing line size, so
            //    the largest-overlap branch is always explored first.
            for (int lid : incidence_map[best_p]) {
                const Line& L  = active_lines[lid];

                // new_bits = uncovered points in this line
                Bitset1024  nb = L.mask.andnot(covered);
                const int npts = nb.count();
                if (npts < 3) continue;

                // Potential loss: sum of point_potential for newly covered pts
                int ploss = 0;
                for (int wi = 0; wi < 16; ++wi) {
                    uint64_t w = nb.words[wi];
                    while (w) {
                        ploss += point_potential[(wi << 6) + std::countr_zero(w)];
                        w &= w - 1;
                    }
                }

                // Apply + recurse + rollback (no 128 B copy)
                covered |= nb;
                (*this)(cur_gain + (npts - 2), cur_pot - ploss);
                covered.andnot_inplace(nb);
            }

            // ── Branch skip: defer best_p to pair coverage ───────────────
            //    Mark covered, recurse, rollback with clear_bit (saves 128 B copy)
            covered.set(best_p);
            (*this)(cur_gain, cur_pot - best_p_pp);
            covered.clear_bit(best_p);
        }
    };

    PointCentricDFS dfs{best_gain, active_pts_vec, incidence_map,
                        active_lines, point_potential};
    dfs(0, total_potential);

    return best_gain;
}

// ═══════════════════════════════════════════════════════════════════════════
// §4  Incremental Solver
//     KEY CHANGE: union-find replaces adj vector + BFS in recompute()
// ═══════════════════════════════════════════════════════════════════════════

class IncrementalSolver {
public:
    static constexpr int MAXN = 1024;

    explicit IncrementalSolver(const std::vector<CollinearLine>& lines)
        : lines_(lines)
        , cur_size_(lines.size(), 0)
    {
        for (int li = 0; li < static_cast<int>(lines.size()); ++li) {
            const auto& L = lines[static_cast<std::size_t>(li)];
            register_event(L.activate_at, li, 3);
            for (int k = 3; k < static_cast<int>(L.pts.size()); ++k)
                register_event(L.pts[static_cast<std::size_t>(k)] + 1, li, k + 1);
        }
    }

    [[nodiscard]] int advance(int n)
    {
        if (n >= 1 && n <= MAXN) {
            for (auto [li, new_sz] : events_[static_cast<std::size_t>(n)]) {
                if (cur_size_[static_cast<std::size_t>(li)] == 0)
                    active_ids_.push_back(li);
                cur_size_[static_cast<std::size_t>(li)] = new_sz;
                dirty_ = true;
            }
        }
        if (!dirty_) return cached_G_star_;
        cached_G_star_ = recompute();
        dirty_ = false;
        return cached_G_star_;
    }

private:
    void register_event(int n, int li, int new_sz)
    {
        if (n >= 1 && n <= MAXN)
            events_[static_cast<std::size_t>(n)].emplace_back(li, new_sz);
    }

    [[nodiscard]] int recompute()
    {
        if (active_ids_.empty()) return 0;

        const int M = static_cast<int>(active_ids_.size());

        // Build per-active-line pts slices and masks
        std::vector<std::vector<int>> lpts (static_cast<std::size_t>(M));
        std::vector<Bitset1024>       lmasks(static_cast<std::size_t>(M));

        for (int i = 0; i < M; ++i) {
            const int  li  = active_ids_[static_cast<std::size_t>(i)];
            const int  sz  = cur_size_[static_cast<std::size_t>(li)];
            const auto& L  = lines_[static_cast<std::size_t>(li)];

            lpts[static_cast<std::size_t>(i)].assign(
                L.pts.begin(), L.pts.begin() + sz);

            Bitset1024& mask = lmasks[static_cast<std::size_t>(i)];
            for (int k = 0; k < sz; ++k)
                mask.set(L.pts[static_cast<std::size_t>(k)]);
        }

        // ── Union-Find component decomposition ────────────────────────────
        // Replaces explicit adj vector + BFS → same O(M²) intersection work
        // but O(M) space, no vector-of-vector allocation.
        std::vector<int> uf(static_cast<std::size_t>(M));
        std::iota(uf.begin(), uf.end(), 0);

        // Path-halving find
        auto find = [&](int x) {
            while (uf[static_cast<std::size_t>(x)] != x) {
                uf[static_cast<std::size_t>(x)] =
                    uf[static_cast<std::size_t>(uf[static_cast<std::size_t>(x)])];
                x = uf[static_cast<std::size_t>(x)];
            }
            return x;
        };

        for (int i = 0; i < M; ++i)
            for (int j = i + 1; j < M; ++j)
                if (lmasks[static_cast<std::size_t>(i)].intersects(
                        lmasks[static_cast<std::size_t>(j)])) {
                    int ri = find(i), rj = find(j);
                    if (ri != rj) uf[static_cast<std::size_t>(ri)] = rj;
                }

        // Gather members per root
        std::unordered_map<int, std::vector<int>> by_root;
        by_root.reserve(static_cast<std::size_t>(M));
        for (int i = 0; i < M; ++i)
            by_root[find(i)].push_back(i);

        std::vector<std::vector<int>> components;
        components.reserve(by_root.size());
        for (auto& [root, members] : by_root)
            components.push_back(std::move(members));

        // Largest component first → best warm-start placement
        std::ranges::sort(components, [](const auto& a, const auto& b) {
            return a.size() > b.size();
        });

        // Solve each component independently and sum.
        // The first (largest) component receives the cached N-1 answer as its
        // initial lower bound, forwarding the N-1 warm-start from §4 into §3.
        int G_star = 0;
        bool first = true;
        for (auto& comp : components) {
            const int lb = first ? cached_G_star_ : 0;
            first = false;
            G_star += solve_max_gain_component(
                std::span<const std::vector<int>>(lpts),
                std::span<const Bitset1024>(lmasks),
                std::move(comp), lb);
        }
        return G_star;
    }

    const std::vector<CollinearLine>&                      lines_;
    std::array<std::vector<std::pair<int,int>>, MAXN + 1> events_;
    std::vector<int>  cur_size_;
    std::vector<int>  active_ids_;
    int               cached_G_star_ = 0;
    bool              dirty_         = false;
};

// ═══════════════════════════════════════════════════════════════════════════
// §5  main
// ═══════════════════════════════════════════════════════════════════════════

int main()
{
    using Clock = std::chrono::steady_clock;
    constexpr int N = 1024;

    const auto t_total_start = Clock::now();

    std::print("Generating primes…\n");
    const std::vector<int> primes = generate_primes(N);

    std::print("Enumerating collinear groups…\n");
    const auto t_enum_start = Clock::now();
    const std::vector<CollinearLine> all_lines = enumerate_collinear_lines(primes);
    const double enum_sec =
        std::chrono::duration<double>(Clock::now() - t_enum_start).count();

    std::print("Found {} collinear groups of size ≥3  (enumeration: {})\n\n",
               all_lines.size(), fmt_dur(enum_sec));

    IncrementalSolver solver(all_lines);

    std::atomic<bool> keep_running{true};
    std::atomic<int>  active_n{0};
    auto current_n_start = Clock::now();

    // Live UI thread: prints current N and elapsed time at 100 ms cadence
    std::thread ui_thread([&]() {
        while (keep_running) {
            if (active_n > 0) {
                auto now = Clock::now();
                double lap = std::chrono::duration<double>(now - current_n_start).count();
                double tot = std::chrono::duration<double>(now - t_total_start).count();
                std::print("\r\033[K[Computing N={:<4} | Current: {} | Total: {}]",
                           active_n.load(), fmt_dur(lap), fmt_dur(tot));
                std::fflush(stdout);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    });

    for (int n = 1; n <= N; ++n) {
        current_n_start = Clock::now();
        active_n = n;

        const int G_star    = solver.advance(n);
        const int min_lines = (n - G_star + 1) / 2;

        const double elapsed =
            std::chrono::duration<double>(Clock::now() - current_n_start).count();

        std::print("\r\033[KN={:<4}  prime={:<8}  lines={:<4}  time={}\n",
                   n, primes[static_cast<std::size_t>(n - 1)],
                   min_lines, fmt_dur(elapsed));
    }

    keep_running = false;
    if (ui_thread.joinable()) ui_thread.join();

    const double total_sec =
        std::chrono::duration<double>(Clock::now() - t_total_start).count();
    std::print("\nTotal wall-clock time: {}\n", fmt_dur(total_sec));
}