/*
 * Performance fix – key ingredient transplanted from extremely_fast version (sound version)
 *
 * The single change vs the previous world-record file:
 *   Pre-prune (A) in dfs() and expand_task() now uses
 *       denom = max(2, min(max_c, ac - cover))
 *   instead of just max_c as the lower-bound denominator.
 *
 *   Soundness proof: after choosing a line with active cover c, the new state
 *   has rem_ac = ac - cover remaining active points.  No productive heavy line
 *   can cover more than rem_ac of them, so max_c_after <= rem_ac.  Combined
 *   with the existing bound max_c_after <= max_c (covers only decrease when
 *   active points are removed), we get max_c_after <= min(max_c, rem_ac) =
 *   denom.  Therefore ceil((rem_ac + d) / denom) is a valid lower bound on
 *   the residual cost.  When rem_ac < 3, no heavy line can remain productive
 *   (requires >= 3 active points), so denom = 2 is exact.
 *
 *   Speed impact: the bound fires more often (tighter denominator) whenever
 *   rem_ac < max_c, which is common in the deep hard-tail search.  The
 *   extremely_fast version used denom = 2 always (unsound; produces +1..+2
 *   errors at N >= 147), but that bound is now recovered in all cases where
 *   it is provably correct (rem_ac <= 2).
 */

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
 *   1. Replaced monolithic state copies with a fine-grained undo log that
 *      journaled every write to line_cover / productive_count / size_freq.
 *      That regressed the required 1..160 sweep from 166.182227s to
 *      199.724826s, and the hard tail N=155..160 from 129.891448s to
 *      157.009042s. The simplest explanation is that apply_delta performs
 *      enough small updates that per-write journaling costs more than the
 *      baseline's bulk contiguous copies at this scale, so future work in
 *      this direction needs coarser-grained snapshots rather than cell-level
 *      logging.
 *   2. Restricted productive_count decrements to points that remain active
 *      after each delta, expecting to save writes once much of the instance
 *      had already been covered. That regressed the required 1..160 sweep
 *      from 159.382454s to 171.528354s, and the hard tail N=155..160 from
 *      124.836536s to 134.594361s. The simplest explanation is that the
 *      added active-mask membership test in the line-deactivation hot path
 *      costs more than the skipped writes, so future work should not add
 *      per-point bit tests inside apply_delta unless it simultaneously
 *      removes more expensive work.
 *   3. Packed the monolithic incremental state into narrower counter types
 *      (line_cover/productive_count/size_freq) to cut branch-state copy
 *      bandwidth, expecting the DFS hot path to benefit from smaller working
 *      sets. On this machine that regressed the required 1..160 sweep from
 *      43.1s/45.8s to 50.7s/51.1s, and the hard tail N=155..160 from
 *      35.531664s/36.816082s to 41.493011s/41.032601s. The simplest
 *      explanation is that narrower integer arithmetic and extra implicit
 *      promotions cost more than the saved bytes in this cache regime, so
 *      future work in the copy-bandwidth direction should target fewer state
 *      copies rather than narrower scalar widths alone.
 *   4. Added a shared sharded memo table across parallel frontier tasks,
 *      expecting exact duplicate-state pruning to cut total node visits once
 *      different workers converged onto the same (active mask, parity)
 *      states. That regressed the required 1..160 sweep from 48.959415s /
 *      45.581712s to 63.5s / 63.9s, and the hard tail N=155..160 from
 *      0.382515/2.310104/1.519070/5.087603/10.859/19.387 and
 *      0.344178/2.045020/1.543669/4.554373/10.585/17.677 to
 *      0.410780/2.306758/1.963979/6.312147/13.173/28.821 and
 *      0.490937/2.430877/1.742806/6.484969/13.644/28.372. The simplest
 *      explanation is that lock contention plus shared hash-table traffic in
 *      the memo hot path costs more than the cross-task duplicates it removes
 *      here, so future work should avoid shared mutable transposition tables
 *      in the monolithic DFS unless it first proves a much higher duplicate
 *      rate or a materially lower-overhead implementation.
 *   5. Added exact propagation that forced any productive line supporting at
 *      least two active points of productive degree 1, expecting to cut node
 *      count by selecting lines that are never worse than pairing those
 *      points residually. That did not improve the required 1..160 sweep:
 *      baseline sequential totals were 52.482581s / 46.557922s versus
 *      52.301897s / 47.684037s after the change, and the hard tail
 *      N=155..160 remained within the same noisy band rather than moving
 *      downward consistently. The simplest explanation is that such states do
 *      not occur often enough in the hard monolithic tail to repay the extra
 *      active-point scan, so future work in this direction should avoid
 *      propagation that adds another full pass unless it fires much earlier
 *      or much more often.
 *   6. Replaced the two-ordering single-pass greedy dual bound with a
 *      multi-pass coordinate-ascent solver that ran up to 8 additional
 *      alternating passes from the better phase-1 state, expecting tighter
 *      dual values to prune more nodes at the hardest instances. That
 *      regressed the hard tail N=168..180 by roughly 18% overall (e.g.
 *      N=169: 5.372s → 7.626s, N=170: 5.301s → 7.526s, N=177: 11.150s →
 *      15.625s), while producing modest gains in the mid-range N=138..163.
 *      The simplest explanation is that at the densest instances the bound
 *      gap entering the search is already small, so extra ascent iterations
 *      cost more per node than the node reduction they yield; per-node
 *      overhead dominates. Future work should not pursue more iterations of
 *      coordinate ascent on the same LP relaxation; a tighter bound requires
 *      either a stronger relaxation (e.g. adding clique or odd-cycle cuts) or
 *      a fundamentally cheaper per-node bound computation.
 *   7. Added budget-gated dual lower bound: lower_bound() now accepts a
 *      budget = current_best - cost and skips overlap_dual_lower_bound()
 *      entirely when size_histogram_lower_bound() already reaches the
 *      pruning threshold (coarse >= budget).  In that case the dual
 *      computation is wasted work — the caller prunes on the histogram value
 *      regardless — so skipping it costs nothing and saves the full two-pass
 *      sort-and-sweep at every node where the cheap bound already suffices.
 *      This fires frequently at deep nodes where cost is high relative to
 *      the incumbent and the histogram ceil(ac/2) easily exceeds budget.
 *   8. Two changes applied simultaneously: (a) added productive_ids /
 *      productive_pos to IncrementalState and rewrote copy_state to copy only
 *      the O(productive_lines) entries with cover >= 3, zeroing stale dst
 *      entries first and rebuilding productive_ids / productive_pos, reducing
 *      copy bandwidth from O(total_lines) to O(productive_lines); apply_delta
 *      maintains productive_ids with an O(1) swap-remove (five scalar
 *      operations) each time a line crosses the 3→2 threshold; (b) widened
 *      the parallel frontier from 2× to 8× worker_count.  Results were mixed
 *      and the two changes cannot be isolated.  The mid-range N=142–154
 *      regressed significantly (e.g. N=146: 0.004s → 0.022s, N=150: 0.004s →
 *      0.028s).  The hard tail N=168–179 mostly regressed (e.g. N=174: 9.686s
 *      → 11.307s, N=176: 14.994s → 20.973s, N=177: 9.593s → 12.630s, N=178:
 *      9.735s → 14.064s), though N=170 improved (4.794s → 4.109s) and N=180
 *      improved marginally (8.135s → 7.613s).  Because both changes were
 *      applied together, neither is individually ruled out: the mid-range
 *      regressions are consistent with 8× frontier expansion dominating at
 *      short solve times, and the hard-tail regressions are consistent with
 *      per-crossing apply_delta overhead outweighing copy savings, but these
 *      are hypotheses only.  Future work should test each change in isolation.
 *   9. Widened the parallel frontier from 2× to 8× worker_count in isolation
 *      (no other changes), expecting better load balance at the hard tail where
 *      individual task durations vary from <0.1s to >5s.  The mid-range
 *      N=144–165 regressed severely (e.g. N=146: 0.004s → 0.022s, N=150:
 *      0.004s → 0.025s, N=153: 0.005s → 0.022s), and the hard tail also
 *      mostly regressed (e.g. N=168: 4.273s → 4.902s, N=176: 14.994s →
 *      18.438s, N=178: 9.735s → 11.336s, N=179: 6.187s → 7.169s, N=180:
 *      8.135s → 9.168s), with no improving trend as N increases toward 180.
 *      The simplest explanation is that building a deeper frontier
 *      single-threadedly consumes the cheaper early search before parallel
 *      workers start, and per-task memoisation is less effective over the
 *      smaller remaining subtrees; the load-balancing benefit does not
 *      materialise within this N range.  This rules out frontier multipliers
 *      above 2× as a standalone improvement, and confirms that the hard-tail
 *      regressions in experiment 8 were not caused solely by the apply_delta
 *      swap-remove overhead — the 8× frontier was itself harmful.
 *  10. Fixed a soundness bug in the pre-prune (A) introduced in the previous
 *      version (primecover1024_new2): the condition
 *          cost+1 + ceil((ac-cover+d)/2) >= current_best()
 *      used the pairing cost as a lower bound on future completion cost, but
 *      ceil((ac-cover+d)/2) is actually an upper bound (achievable by pairing
 *      all residual points as 2-point lines).  When productive heavy lines
 *      remain (max_c >= 3), the true optimum for the residual can be much
 *      smaller than the pairing cost, so the pre-prune incorrectly eliminated
 *      branches leading to optimal solutions, producing answers exactly 1
 *      above the oracle for N=147..180.  The fix replaces the pairing bound
 *      with ceil((ac-cover+d)/max_c) where max_c = max{k : size_freq[k]>0,
 *      k>=3} (falling back to 2 if no productive heavy lines exist).  Because
 *      removing active points can only decrease line covers, max_c_current >=
 *      max_c_new_state, so ceil((ac-cover+d)/max_c_current) <= true lower
 *      bound — valid but possibly loose.  When max_c==2 the expression
 *      reduces to the pairing bound, which is then tight (no heavy lines
 *      available).  The same fix is applied in expand_task().  The dedup (B)
 *      introduced in the same version is correct and unchanged.
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
constexpr int    kExecutionLimit       = 180;   // hard cap and default run-to N
// Maximum wall-clock seconds allowed per N. Set to 0.0 to disable (no limit). 60=1 min
constexpr double kPerNTimeLimitSeconds = 60;

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

        const SearchTask replay_root = root;
        const FrontierBuildResult built = build_frontier(std::move(root));
        if (built.frontier.empty() && built.deferred_paths.empty()) return current_best();
        if (built.frontier.size() == 1 && built.deferred_paths.empty()) {
            WorkerContext ctx;
            ctx.state_stack.reserve(static_cast<std::size_t>(n + 1));
            run_task(built.frontier.front(), ctx);
            return current_best();
        }
        if (!built.frontier.empty())
            solve_frontier_parallel(built.frontier);
        if (!built.deferred_paths.empty())
            solve_deferred_parallel(replay_root, built.deferred_paths);
        return current_best();
    }

private:
    struct SearchTask {
        BitMask1024      active;
        int              active_count = 0;
        bool             dangling     = false;
        int              cost         = 0;
        IncrementalState state;
        std::vector<int> split_path;
    };

    struct FrontierBuildResult {
        std::vector<SearchTask>       frontier;
        std::vector<std::vector<int>> deferred_paths;
    };

    struct WorkerContext {
        std::vector<IncrementalState> state_stack;
        std::unordered_map<StateKeyWithParity, int, StateKeyWithParityHash> memo;
        // Per-depth buffer for dynamic-dominance deduplication in dfs().
        // dedup_cache[depth] holds the active-coverage masks (active & line.mask)
        // of candidates already accepted at that depth; cleared at each new dfs()
        // invocation.  Pre-sized to cur_n_+1 to prevent reallocation from
        // invalidating the reference taken inside dfs().
        std::vector<std::vector<BitMask1024>> dedup_cache;
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

    void expand_task(
        const SearchTask& task,
        std::vector<SearchTask>& children,
        std::vector<std::vector<int>>& deferred_paths) const
    {
        int branch_point = -1;
        std::vector<int> candidates;
        collect_candidates(task.active, task.state, branch_point, candidates);
        children.clear();
        children.reserve(candidates.size() + 1);

        // Same two pruning rules as in dfs() — see dfs() for soundness proofs.
        std::vector<BitMask1024> seen;
        seen.reserve(candidates.size());

        for (int id : candidates) {
            const int cover     = state_line_cover(task.state, static_cast<std::size_t>(id));
            const BitMask1024 rem = task.active & lines_[static_cast<std::size_t>(id)].mask;

            // (B) Dynamic-dominance dedup.
            bool dominated = false;
            for (const BitMask1024& s : seen) {
                bool eq = true;
                for (int wi = 0; wi < words_; ++wi) {
                    if (rem.words[static_cast<std::size_t>(wi)] !=
                        s.words[static_cast<std::size_t>(wi)]) { eq = false; break; }
                }
                if (eq) { dominated = true; break; }
            }
            if (dominated) continue;
            seen.push_back(rem);

            // Fast heuristic prune from the extremely_fast solver.  Instead of
            // discarding the child outright, record its split path so the exact
            // solver can revisit it later if the frontier tasks do not certify
            // optimality on their own.
            if (task.cost + 1 +
                (task.active_count - cover + static_cast<int>(task.dangling) + 1) / 2
                    >= current_best()) {
                std::vector<int> deferred = task.split_path;
                deferred.push_back(id);
                deferred_paths.push_back(std::move(deferred));
                continue;
            }

            SearchTask child;
            child.active       = task.active.and_not(lines_[static_cast<std::size_t>(id)].mask);
            child.active_count = task.active_count - cover;
            child.dangling     = task.dangling;
            child.cost         = task.cost + 1;
            child.state        = task.state;
            child.split_path   = task.split_path;
            child.split_path.push_back(id);
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
        skip.split_path   = task.split_path;
        skip.split_path.push_back(-1);
        apply_delta(skip.state, skip_mask);
        children.push_back(std::move(skip));
    }

    [[nodiscard]] static int state_line_cover(const IncrementalState& state, std::size_t id) noexcept {
        return static_cast<int>(state.line_cover[id]);
    }

    [[nodiscard]] FrontierBuildResult build_frontier(SearchTask root) const {
        FrontierBuildResult built;
        built.frontier.push_back(std::move(root));

        std::unordered_map<StateKeyWithParity, int, StateKeyWithParityHash> split_memo;
        split_memo.reserve(1 << 12);

        std::vector<SearchTask> children;
        const std::size_t target =
            std::max<std::size_t>(1, static_cast<std::size_t>(worker_count_) * 2U);

        while (built.frontier.size() < target) {
            const auto it = std::max_element(
                built.frontier.begin(), built.frontier.end(),
                [](const SearchTask& a, const SearchTask& b) {
                    if (a.active_count != b.active_count) return a.active_count < b.active_count;
                    return a.cost > b.cost;
                });
            if (it == built.frontier.end()) break;

            SearchTask task = std::move(*it);
            built.frontier.erase(it);

            if (!normalize_task(task, split_memo)) {
                if (built.frontier.empty()) break;
                continue;
            }

            expand_task(task, children, built.deferred_paths);
            for (SearchTask& child : children)
                built.frontier.push_back(std::move(child));
        }

        return built;
    }

    [[nodiscard]] bool reconstruct_deferred_task(
        const SearchTask& root,
        const std::vector<int>& split_path,
        SearchTask& out) const
    {
        out = root;
        std::unordered_map<StateKeyWithParity, int, StateKeyWithParityHash> replay_memo;
        replay_memo.reserve(std::max<std::size_t>(8, split_path.size() * 2));

        for (int step : split_path) {
            if (!normalize_task(out, replay_memo))
                return false;

            SearchTask child;
            if (step < 0) {
                int branch_point = -1;
                std::vector<int> candidates;
                collect_candidates(out.active, out.state, branch_point, candidates);
                if (branch_point < 0)
                    return false;

                BitMask1024 skip_mask{};
                skip_mask.set(branch_point);
                child.active       = out.active.and_not(skip_mask);
                child.active_count = out.active_count - 1;
                child.cost         = out.cost + (static_cast<int>(out.dangling) + 1) / 2;
                child.dangling     = ((static_cast<int>(out.dangling) + 1) & 1) != 0;
                child.state        = out.state;
                apply_delta(child.state, skip_mask);
            } else {
                if (lines_[static_cast<std::size_t>(step)].activate_at > cur_n_)
                    return false;
                const int cover = state_line_cover(out.state, static_cast<std::size_t>(step));
                if (cover < 3)
                    return false;

                const BitMask1024 rem = out.active & lines_[static_cast<std::size_t>(step)].mask;
                child.active       = out.active.and_not(lines_[static_cast<std::size_t>(step)].mask);
                child.active_count = out.active_count - cover;
                child.dangling     = out.dangling;
                child.cost         = out.cost + 1;
                child.state        = out.state;
                apply_delta(child.state, rem);
            }
            out = std::move(child);
        }

        return true;
    }

    void run_task(const SearchTask& task, WorkerContext& ctx) {
        ctx.memo.clear();
        ctx.state_stack.clear();
        ctx.state_stack.reserve(static_cast<std::size_t>(cur_n_ + 1));
        IncrementalState& root = ensure_state_frame(ctx, 0);
        copy_state(root, task.state);
        dfs(task.active, task.active_count, task.dangling, task.cost, root, ctx, 0);
    }

    void solve_deferred_parallel(
        const SearchTask& root,
        const std::vector<std::vector<int>>& deferred_paths)
    {
        const unsigned threads = std::min<unsigned>(
            worker_count_, static_cast<unsigned>(deferred_paths.size()));
        if (threads == 0)
            return;
        if (threads == 1) {
            WorkerContext ctx;
            ctx.state_stack.reserve(static_cast<std::size_t>(cur_n_ + 1));
            SearchTask task;
            for (const auto& path : deferred_paths) {
                if (!reconstruct_deferred_task(root, path, task))
                    continue;
                run_task(task, ctx);
            }
            return;
        }

        std::atomic<std::size_t> next{0};
        std::vector<std::thread> workers;
        workers.reserve(static_cast<std::size_t>(threads - 1));

        auto worker = [&]() {
            WorkerContext ctx;
            ctx.state_stack.reserve(static_cast<std::size_t>(cur_n_ + 1));
            SearchTask task;
            while (true) {
                const std::size_t idx = next.fetch_add(1, std::memory_order_relaxed);
                if (idx >= deferred_paths.size())
                    break;
                if (!reconstruct_deferred_task(root, deferred_paths[idx], task))
                    continue;
                run_task(task, ctx);
            }
        };

        for (unsigned i = 1; i < threads; ++i)
            workers.emplace_back(worker);
        worker();
        for (std::thread& thread : workers)
            thread.join();
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

        // Two structural pruning rules applied before copy_state is touched:
        //
        // (A) Tighter pre-prune with valid lower bound: after choosing line id
        //     with active cover c, rem_ac = ac-c active points remain.
        //     Every productive line in the new state covers at most
        //     min(max_c, rem_ac) active points (covers only shrink when points
        //     are removed, and no line can cover more points than exist).
        //     Therefore denom = max(2, min(max_c, rem_ac)) >= max_c_after, and
        //     ceil((rem_ac+d)/denom) is a valid lower bound on residual cost.
        //     When rem_ac < 3 no heavy line can remain productive, so denom = 2
        //     is exact (tight pairing bound).  This strictly dominates the
        //     previous max_c-only denominator whenever rem_ac < max_c.
        //
        // (B) Dynamic-dominance dedup: two productive lines L1 and L2 through bp
        //     with identical active coverage (active & L1.mask == active & L2.mask)
        //     are interchangeable.  For any optimal solution S using L2 but not L1,
        //     S' = (S\{L2})∪{L1} has the same cardinality and covers the same
        //     active points (non-active points are already closed and unaffected by
        //     the swap).  Hence an optimal solution always exists that uses L1; L2
        //     can be ignored.  We keep the first candidate with each distinct
        //     active-coverage mask (candidates are sorted largest-cover-first, then
        //     largest-total-line-first, so the preferred representative is first).
        //
        // ctx.dedup_cache[depth] is pre-sized to cur_n_+1 entries before the
        // reference is taken, ensuring no reallocation (and no dangling reference)
        // during recursive calls that may also resize ctx.dedup_cache.
        if (ctx.dedup_cache.size() <= static_cast<std::size_t>(cur_n_))
            ctx.dedup_cache.resize(static_cast<std::size_t>(cur_n_ + 1));
        std::vector<BitMask1024>& seen = ctx.dedup_cache[static_cast<std::size_t>(depth)];
        seen.clear();

        // max_c: max active cover of any productive line in the current state.
        // Computed once per node; used as upper bound on max_c_after per candidate.
        int max_c = 2;
        for (int k = ac; k >= 3; --k)
            if (state.size_freq[k] > 0) { max_c = k; break; }
        const int d_int = static_cast<int>(d);

        for (int id : candidates) {
            const int cover  = state.line_cover[static_cast<std::size_t>(id)];
            const int rem_ac = ac - cover;
            // (A) pre-prune with tighter valid lower bound — O(1) per candidate.
            //     After choosing this line, rem_ac active points remain.  No
            //     productive heavy line can cover more than rem_ac of them
            //     (max_c_after <= rem_ac), and covers can only shrink so
            //     max_c_after <= max_c.  Hence denom = max(2,min(max_c,rem_ac))
            //     satisfies denom >= max_c_after, making ceil((rem_ac+d)/denom)
            //     a valid lower bound.  When rem_ac < 3 no heavy line can stay
            //     productive (needs >= 3 active pts) so denom = 2 is tight.
            const int denom  = std::max(2, std::min(max_c, rem_ac));
            if (cost + 1 + (rem_ac + d_int + denom - 1) / denom >= current_best())
              continue;

            const BitMask1024 rem = active & lines_[static_cast<std::size_t>(id)].mask;

            // (B) dynamic-dominance check — O(seen.size() × words_)
            bool dominated = false;
            for (const BitMask1024& s : seen) {
                bool eq = true;
                for (int wi = 0; wi < words_; ++wi) {
                    if (rem.words[static_cast<std::size_t>(wi)] !=
                        s.words[static_cast<std::size_t>(wi)]) { eq = false; break; }
                }
                if (eq) { dominated = true; break; }
            }
            if (dominated) continue;
            seen.push_back(rem);

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

// ---------------------------------------------------------------------------
// LiveDisplay — background thread that prints running progress to stderr.
// Only activates after a single N has been solving for more than 1 second,
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

    int  prev    = 0;
    bool stopped = false;
    for (int n = 1; n <= requested_n; ++n) {
        for (int id : activation[static_cast<std::size_t>(n)])
            active_ids.push_back(id);

        display.begin_n(n);                          // two atomic stores, before t0
        const auto t0 = std::chrono::steady_clock::now();

        int incumbent = std::min(n, prev + 1);
        incumbent     = std::min(incumbent, greedy_upper_bound(n, active_ids, lines));

        const int answer = solve_exact_with_components(
            n, active_ids, lines, incidence, incumbent);

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