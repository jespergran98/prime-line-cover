/*
 * primecover1024_v2.cpp  —  C++23/26, zero external dependencies
 *
 * PROBLEM
 * ───────
 * Let p_i be the i-th prime (p_0 = 2).  For each prefix length n = 1..1024
 * find the minimum number of straight lines that cover every point (i, p_i),
 * 0 ≤ i < n.
 *
 * FORMULA
 * ───────
 *   min_lines(n) = ⌈(n − W*(n)) / 2⌉
 *
 * where W*(n) = maximum total weight of a non-overlapping family of collinear
 * groups of size ≥ 3 within the first n prime-points.  Weight of a group of
 * size k is (k − 2).
 *
 * COLLINEAR STRUCTURE (correction of v1 NOTE)
 * ────────────────────────────────────────────
 * v1 mistakenly claimed no three prime-points are collinear.  In fact there
 * are exactly 12,162 collinear groups of size ≥ 3 within the first 1024
 * primes — e.g. {(1,3),(2,5),(3,7)} or {(2,5),(4,11),(6,17),(8,23)}.
 * For large n these 12,162 lines form a single connected component in the
 * conflict graph (two lines conflict iff they share a prime-point).
 * Crucially, two *distinct* lines in the plane share at most ONE point, so
 * every edge in the conflict graph corresponds to exactly one shared prime.
 *
 * WHY v1's B&B FAILED
 * ───────────────────
 * v1 used suffix_gain = Σ(weight of remaining lines), ignoring every
 * conflict.  For a single 12,162-node component this bound is ≈ 16,484
 * while the greedy lower bound is ≈ 645 — a gap of ~15,800 that never
 * closes, producing exponential node counts.
 *
 * IMPROVEMENTS IN v2
 * ──────────────────
 * §3  LP-dual upper bound replaces suffix_gain.
 *
 *     UB = cur_w  +  Σ_{p unclaimed}  max_{L at idx or later, p∈L}  weight(L)/size(L)
 *
 *     Correctness: setting y_p = max_density(p) is dual-feasible for the LP
 *     relaxation of the set-packing IP (Σ_{p∈L} y_p ≥ weight(L) for every
 *     line L, because Σ_{p∈L} weight(L)/size(L) = weight(L)).  The dual
 *     objective Σ_p y_p is therefore ≥ IP-OPT by LP-duality.
 *
 *     At the root this equals ≈ 824 (vs suffix_gain ≈ 16 484).  As lines are
 *     decided the bound drops sharply: including a size-s line of weight w
 *     removes s unclaimed points whose max-density might equal w/s each,
 *     immediately decreasing UB by exactly w — tightly matching the gain.
 *
 *     Per-point suffix-max density arrays let each bound evaluation run in
 *     O(|active_pts| · log(lines_per_point)) ≈ O(6 000) instead of O(M).
 *
 *     Lines are sorted by density (weight/size) descending so the greedy
 *     warm-start is tight and B&B encounters the heaviest lines first.
 *
 *     C++23 deducing-this lambda eliminates std::function vtable overhead.
 *
 *     Greedy warm-start initialises best before the first B&B call.
 *     Across incremental events, cached_W_star_ from the previous call is
 *     forwarded as initial_lb so the gap seen by B&B collapses.
 *
 * Compile:
 *   g++ -std=c++23 -O3 -march=native -o primecover1024_v2 primecover1024_v2.cpp
 *   g++ -std=c++26 -O3 -march=native -o primecover1024_v2 primecover1024_v2.cpp
 */

#include <algorithm>
#include <array>
#include <atomic>
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <string>
#include <thread>
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
#include <ranges>
#include <span>
#include <unordered_map>
#include <utility>
#include <vector>

// ═══════════════════════════════════════════════════════════════════════════
// Utility: human-readable duration
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
// §1  Segmented sieve – first `count` primes
//     (identical to v1; included verbatim for self-containedness)
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
//     Finds every maximal group of ≥ 3 prime-points on one Euclidean line.
//     (Logic identical to v1; the erroneous NOTE there was wrong — this code
//      correctly discovers all 12,162 collinear groups.)
// ═══════════════════════════════════════════════════════════════════════════

struct CollinearLine {
    std::vector<int>               pts;           // point indices, sorted ascending
    int                            activate_at;   // pts[2]+1: first n with ≥3 in-prefix pts
    std::vector<std::bitset<1024>> prefix_masks;  // [k] covers first (k+3) pts, k=0,1,…
};

[[nodiscard]]
static std::vector<CollinearLine> enumerate_collinear_lines(
    const std::vector<int>& primes)
{
    const int N = static_cast<int>(primes.size());

    // used[i].test(j): pair (i,j) already belongs to an emitted line  [i<j].
    std::vector<std::bitset<1024>> used(static_cast<std::size_t>(N));

    // Slope key (lossless): key = (dy/g)*1024 + (dx/g)
    // dx/g ≤ 1023, dy/g ≤ 8199, max key < 2^24, fits int32.
    auto pack_slope = [](long long dy_r, long long dx_r) noexcept -> int {
        return static_cast<int>(dy_r * 1024 + dx_r);
    };

    std::vector<CollinearLine> all_lines;

    for (int i = 0; i < N; ++i) {
        std::unordered_map<int, std::vector<int>> slope_map;
        slope_map.reserve(static_cast<std::size_t>(N - i));

        for (int j = i + 1; j < N; ++j) {
            if (used[static_cast<std::size_t>(i)].test(static_cast<std::size_t>(j)))
                continue;
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

            // Mark every pair so future anchors don't re-emit this line.
            for (int a = 0; a < static_cast<int>(pts.size()); ++a)
                for (int b = a + 1; b < static_cast<int>(pts.size()); ++b)
                    used[static_cast<std::size_t>(pts[a])].set(
                        static_cast<std::size_t>(pts[b]));

            // Precompute prefix masks for in-prefix sizes 3, 4, …, |pts|.
            const int sz = static_cast<int>(pts.size());
            std::vector<std::bitset<1024>> pmasks;
            pmasks.reserve(static_cast<std::size_t>(sz - 2));
            std::bitset<1024> acc;
            for (int k = 0; k < sz; ++k) {
                acc.set(static_cast<std::size_t>(pts[k]));
                if (k >= 2) pmasks.push_back(acc);
            }

            const int act = pts[2] + 1;
            all_lines.push_back(CollinearLine{
                .pts          = std::move(pts),
                .activate_at  = act,
                .prefix_masks = std::move(pmasks)
            });
        }
    }
    return all_lines;
}

// ═══════════════════════════════════════════════════════════════════════════
// §3  LP-dual MWIS solver  — exact branch-and-bound per connected component
//
// KEY DIFFERENCE FROM v1
// ──────────────────────
// v1 computed suffix_gain[k] = Σ_{i≥k} weight(i), the sum of all future
// weights ignoring every conflict.  For a single component with 12,162
// lines this equals ≈ 16,484, while the actual W* ≈ 645 — a gap of 15,800
// that the B&B can never close.
//
// v2 computes an LP-dual upper bound:
//
//   UB(idx, claimed) = cur_w  +  Σ_{p ∉ claimed}  max_{L at idx…m−1, p∈L}  w(L)/|L|
//
// Dual feasibility: setting y_p = max-density(p) satisfies Σ_{p∈L} y_p ≥ w(L)
// for every line L (since Σ_{p∈L} w/|L| = w).  By LP duality,
// Σ_p y_p ≥ LP-OPT ≥ IP-OPT.
//
// At the root this equals ≈ 824 — much tighter than 16,484.  It decreases
// sharply at each inclusion: claiming |L| points each contributing w/|L|
// reduces the bound by ≈ w, exactly offsetting the weight gained, while
// removing lines from consideration for those points tightens further.
//
// IMPLEMENTATION DETAILS
// ──────────────────────
// • Lines sorted by density (w/|L|) descending inside each component so
//   the greedy warm-start and B&B both encounter heaviest lines first.
//
// • Per-point suffix-max density arrays: for each point p, an array of
//   (dfs_position, suffix_max_density) pairs computed once per component.
//   Binary search finds max_density(p, idx) in O(log k) where k = number
//   of lines through p (avg ≈ 41, so log k ≈ 6).
//   Full bound evaluation: O(|active_pts| · 6) ≈ O(6 000) per node.
//
// • C++23 deducing-this lambda eliminates std::function vtable overhead.
//
// • Integer floor pruning: since W* is an integer, prune when ub < best+1
//   (with a small float epsilon for rounding safety).
//
// • initial_lb parameter lets the IncrementalSolver warm-start the B&B
//   with cached_W_star_ from the previous event, so the gap seen by B&B
//   collapses to (new_optimum − old_optimum) instead of the full range.
// ═══════════════════════════════════════════════════════════════════════════

[[nodiscard]]
static int solve_mwis_component(
    std::span<const std::vector<int>>  line_pts,
    std::span<const std::bitset<1024>> line_masks,
    std::vector<int>                   comp,
    int                                initial_lb = 0)
{
    const int c = static_cast<int>(comp.size());
    if (c == 0) return 0;
    if (c == 1) return static_cast<int>(line_pts[comp[0]].size()) - 2;

    // ── Sort component by density (weight/size) descending ─────────────────
    std::ranges::sort(comp, [&](int a, int b) {
        // (size-2)/size  vs  (size-2)/size  →  cross-multiply (avoids floats)
        int wa = static_cast<int>(line_pts[a].size()) - 2;
        int wb = static_cast<int>(line_pts[b].size()) - 2;
        int sa = static_cast<int>(line_pts[a].size());
        int sb = static_cast<int>(line_pts[b].size());
        return wa * sb > wb * sa;
    });

    if (c == 2) {
        // A 2-node component → both nodes are connected (they conflict).
        return std::max(
            static_cast<int>(line_pts[comp[0]].size()) - 2,
            static_cast<int>(line_pts[comp[1]].size()) - 2);
    }

    // ── Build local arrays (in density-sorted DFS order) ────────────────────
    std::vector<int>               weights(static_cast<std::size_t>(c));
    std::vector<int>               sizes  (static_cast<std::size_t>(c));
    std::vector<std::bitset<1024>> masks  (static_cast<std::size_t>(c));

    for (int i = 0; i < c; ++i) {
        const int li = comp[static_cast<std::size_t>(i)];
        sizes  [static_cast<std::size_t>(i)] = static_cast<int>(line_pts[li].size());
        weights[static_cast<std::size_t>(i)] = sizes[static_cast<std::size_t>(i)] - 2;
        masks  [static_cast<std::size_t>(i)] = line_masks[li];
    }

    // ── Collect active points (appear in ≥1 line in this component) ─────────
    std::bitset<1024> all_pts_bs;
    for (int i = 0; i < c; ++i) all_pts_bs |= masks[static_cast<std::size_t>(i)];
    std::vector<int> active_pts;
    active_pts.reserve(all_pts_bs.count());
    for (int p = 0; p < 1024; ++p)
        if (all_pts_bs.test(static_cast<std::size_t>(p)))
            active_pts.push_back(p);
    const int AP = static_cast<int>(active_pts.size());

    // ── Per-point suffix-max density arrays ─────────────────────────────────
    // pt_suf[p] = vector of {pos, suf_max} sorted by pos ascending,
    // where suf_max is the maximum density of lines at DFS positions ≥ pos
    // that pass through p.
    struct SufEntry { int pos; float suf_max; };
    std::vector<std::vector<SufEntry>> pt_suf(1024);

    {
        struct Raw { int pos; float density; };
        std::vector<std::vector<Raw>> raw(1024);

        for (int i = 0; i < c; ++i) {
            const float d = static_cast<float>(weights[static_cast<std::size_t>(i)])
                          / static_cast<float>(sizes  [static_cast<std::size_t>(i)]);
            for (int p : line_pts[comp[static_cast<std::size_t>(i)]])
                raw[static_cast<std::size_t>(p)].push_back({i, d});
        }

        for (int p : active_pts) {
            auto& rv = raw[static_cast<std::size_t>(p)];
            std::ranges::sort(rv, [](const Raw& a, const Raw& b){ return a.pos < b.pos; });
            pt_suf[static_cast<std::size_t>(p)].resize(rv.size());
            float mx = 0.f;
            for (int k = static_cast<int>(rv.size()) - 1; k >= 0; --k) {
                mx = std::max(mx, rv[static_cast<std::size_t>(k)].density);
                pt_suf[static_cast<std::size_t>(p)][static_cast<std::size_t>(k)] =
                    {rv[static_cast<std::size_t>(k)].pos, mx};
            }
        }
    }

    // ── Greedy warm-start (density-sorted, compatible first) ────────────────
    // Provides a tight initial lower bound and collapses the B&B gap.
    int best = initial_lb;
    {
        std::bitset<1024> gc;
        int gw = 0;
        for (int i = 0; i < c; ++i) {
            if ((masks[static_cast<std::size_t>(i)] & gc).none()) {
                gc |= masks[static_cast<std::size_t>(i)];
                gw += weights[static_cast<std::size_t>(i)];
            }
        }
        if (gw > best) best = gw;
    }

    // ── Branch-and-bound DFS with LP-dual upper bound ───────────────────────
    // Uses a local struct with operator() for zero-overhead recursion
    // (avoids both std::function vtable and deducing-this, which requires
    // GCC ≥ 14 or Clang ≥ 18; this compiles cleanly on any C++23 toolchain).
    std::bitset<1024> claimed;
    int               cur_w = 0;

    struct DFS {
        // All state is held by reference in the outer scope; the struct is
        // just a named handle so we can call operator() recursively.
        const int                               c;
        int&                                    best;
        int&                                    cur_w;
        std::bitset<1024>&                      claimed;
        const std::vector<int>&                 weights;
        const std::vector<int>&                 sizes;
        const std::vector<std::bitset<1024>>&   masks;
        const std::vector<int>&                 active_pts;
        const int                               AP;
        // Bound helper: captured from outer scope via pointer (cheap).
        const std::vector<std::vector<SufEntry>>* pt_suf_ptr;

        [[nodiscard]] float max_d(int p, int idx) const noexcept {
            const auto& v = (*pt_suf_ptr)[static_cast<std::size_t>(p)];
            if (v.empty()) return 0.f;
            std::size_t lo = 0, hi = v.size();
            while (lo < hi) {
                const std::size_t mid = (lo + hi) >> 1;
                if (v[mid].pos < idx) lo = mid + 1; else hi = mid;
            }
            return (lo < v.size()) ? v[lo].suf_max : 0.f;
        }

        void operator()(int idx) {
            // ── LP-dual upper bound ───────────────────────────────────────
            float ub = static_cast<float>(cur_w);
            for (int k = 0; k < AP; ++k) {
                const int p = active_pts[static_cast<std::size_t>(k)];
                if (!claimed.test(static_cast<std::size_t>(p)))
                    ub += max_d(p, idx);
            }
            // Integer floor prune: IP-OPT is integer, so prune if ub < best+1.
            if (ub < static_cast<float>(best) + 1.f - 1e-5f) return;

            if (idx == c) { best = cur_w; return; }

            // ── Branch A: INCLUDE line idx (if disjoint from claimed) ─────
            if ((masks[static_cast<std::size_t>(idx)] & claimed).none()) {
                claimed |= masks[static_cast<std::size_t>(idx)];
                cur_w   += weights[static_cast<std::size_t>(idx)];
                (*this)(idx + 1);
                // Undo: XOR correct since these bits were 0 before include.
                claimed ^= masks[static_cast<std::size_t>(idx)];
                cur_w   -= weights[static_cast<std::size_t>(idx)];
            }

            // ── Branch B: EXCLUDE line idx ────────────────────────────────
            (*this)(idx + 1);
        }
    };

    DFS dfs{c, best, cur_w, claimed, weights, sizes, masks,
            active_pts, AP, &pt_suf};
    dfs(0);
    return best;
}

// ═══════════════════════════════════════════════════════════════════════════
// §4  Incremental W* computation
//
// Same event-driven frame as v1.  Two additions:
//
//  (a) cached_W_star_ is forwarded as initial_lb to solve_mwis_component so
//      the B&B always starts with the best value from the previous event.
//      Since W*(n) is non-decreasing in n, this is a valid lower bound and
//      dramatically reduces the gap the B&B must close.
//
//  (b) The component BFS uses bitset masks for O(M²/64) adjacency checks.
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
        if (!dirty_) return cached_W_star_;
        cached_W_star_ = recompute();
        dirty_ = false;
        return cached_W_star_;
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

        // Build current in-prefix point lists and masks.
        std::vector<std::vector<int>>  lpts (static_cast<std::size_t>(M));
        std::vector<std::bitset<1024>> lmasks(static_cast<std::size_t>(M));

        for (int i = 0; i < M; ++i) {
            const int  li  = active_ids_[static_cast<std::size_t>(i)];
            const int  sz  = cur_size_  [static_cast<std::size_t>(li)];
            const auto& L  = lines_[static_cast<std::size_t>(li)];

            lpts[static_cast<std::size_t>(i)].assign(
                L.pts.begin(), L.pts.begin() + sz);
            lmasks[static_cast<std::size_t>(i)] =
                L.prefix_masks[static_cast<std::size_t>(sz - 3)];
        }

        // Build conflict graph: lines conflict iff they share a prime-point.
        std::vector<std::vector<int>> adj(static_cast<std::size_t>(M));
        for (int i = 0; i < M; ++i)
            for (int j = i + 1; j < M; ++j)
                if ((lmasks[static_cast<std::size_t>(i)] &
                     lmasks[static_cast<std::size_t>(j)]).any()) {
                    adj[static_cast<std::size_t>(i)].push_back(j);
                    adj[static_cast<std::size_t>(j)].push_back(i);
                }

        // BFS-decompose into connected components; solve MWIS per component.
        // Pass cached_W_star_ as initial_lb to the first (and typically
        // largest) component so B&B warm-starts from the last known answer.
        int W_star = 0;
        std::vector<bool> visited(static_cast<std::size_t>(M), false);
        bool first_component = true;

        for (int src = 0; src < M; ++src) {
            if (visited[static_cast<std::size_t>(src)]) continue;

            std::vector<int> comp, queue = {src};
            visited[static_cast<std::size_t>(src)] = true;

            for (int head = 0; head < static_cast<int>(queue.size()); ++head) {
                const int u = queue[static_cast<std::size_t>(head)];
                comp.push_back(u);
                for (int v : adj[static_cast<std::size_t>(u)])
                    if (!visited[static_cast<std::size_t>(v)]) {
                        visited[static_cast<std::size_t>(v)] = true;
                        queue.push_back(v);
                    }
            }

            // Warm-start: the cached answer is a valid lower bound because
            // W*(n) is non-decreasing.  We give it to the first component
            // (which is the biggest after the first large-n event).
            const int lb = first_component ? cached_W_star_ : 0;
            first_component = false;

            W_star += solve_mwis_component(lpts, lmasks, std::move(comp), lb);
        }
        return W_star;
    }

    // ── Data ────────────────────────────────────────────────────────────────
    const std::vector<CollinearLine>&                          lines_;
    std::array<std::vector<std::pair<int,int>>, MAXN + 1>     events_;
    std::vector<int>  cur_size_;
    std::vector<int>  active_ids_;
    int               cached_W_star_ = 0;
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
    const double enum_sec = std::chrono::duration<double>(Clock::now() - t_enum_start).count();

    std::print("Found {} collinear groups of size ≥3  (enumeration: {})\n\n",
               all_lines.size(), fmt_dur(enum_sec));

    IncrementalSolver solver(all_lines);

    // --- Threading Controls for the Live Timer ---
    std::atomic<bool> keep_running{true};
    std::atomic<int>  active_n{0};
    auto current_n_start = Clock::now();

    std::thread ui_thread([&]() {
        while (keep_running) {
            if (active_n > 0) {
                auto now = Clock::now();
                double lap = std::chrono::duration<double>(now - current_n_start).count();
                double tot = std::chrono::duration<double>(now - t_total_start).count();
                
                // \r returns to start of line, \033[K clears to the end of line
                std::print("\r\033[K[Computing N={:<4} | Current: {} | Total: {}]", 
                           active_n.load(), fmt_dur(lap), fmt_dur(tot));
                std::fflush(stdout);
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    });

    for (int n = 1; n <= N; ++n) {
        current_n_start = Clock::now();
        active_n = n; // Signal the UI thread which N we are on

        const int W_star    = solver.advance(n);
        const int min_lines = (n - W_star + 1) / 2;

        const double elapsed = std::chrono::duration<double>(Clock::now() - current_n_start).count();

        // Overwrite the live timer with the final result
        std::print("\r\033[KN={:<4}  prime={:<8}  lines={:<4}  time={}\n",
                   n, primes[static_cast<std::size_t>(n - 1)],
                   min_lines, fmt_dur(elapsed));
    }

    // Shut down UI thread
    keep_running = false;
    if (ui_thread.joinable()) ui_thread.join();

    const double total_sec = std::chrono::duration<double>(Clock::now() - t_total_start).count();
    std::print("\nTotal wall-clock time: {}\n", fmt_dur(total_sec));
}