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
 *   3. Reusing an incremental parent dual inside lagrangian_upper_bound() was
 *      supposed to stop rebuilding the bound from scratch by updating only the
 *      touched points, touched productive lines, and their undo state along
 *      the DFS path. Two direct fast-path versions ran faster but changed the
 *      certified answers, and the guarded exact version still took 86.836s to
 *      reach N=730 against the 62.409s baseline to the same cutoff. The
 *      simplest explanation is that the carried dual state lost too much
 *      pruning power once it was no longer fully rebuilt, so future attempts
 *      should avoid parent-state reuse unless the refreshed child bound and
 *      every derived branch cap can be kept both exact and nearly as tight as
 *      the rebuilt bound.
 *   4. Reusing an incremental child bound state was supposed to speed up the
 *      branch-and-bound search by carrying the active-point view and per-line
 *      dual data across include/exclude steps instead of rebuilding them at
 *      each child. Three exact 750 runs with increasingly selective variants
 *      all kept the same answers but still lost to the 157.028s / 156.433s
 *      baseline, finishing at 250.782s, 171.815s, and 174.168s while the hard
 *      tail around N=731-N=742 slowed down. The simplest explanation is that
 *      the extra bookkeeping and refresh work for the carried bound state cost
 *      more than the partial rebuild they saved, so future attempts should
 *      avoid path-level dual-state maintenance unless the full bound and
 *      restore path can be updated with much lower constant overhead.
 *   5. A same-cost local repair step around the previous exact heavy-line set
 *      was supposed to keep N->N+1 from falling back to a fresh global search
 *      by fixing the cover around new-point lines first. Three exact 750
 *      variants all kept the certified answers, but the broad local-pool run
 *      still lost at 173.722s, the larger-pool follow-up landed at 141.320s
 *      while leaving the hard N=742 tail slower than baseline, and the
 *      anchored new-point version fell back to 235.317s against the
 *      109.739s / 170.727s baseline pair. The simplest explanation is that the
 *      pre-search repair either misses the real same-cost fixes or burns too
 *      much exact-search work before the global fallback, so future attempts
 *      should avoid adding another repair layer unless it can prove a carry-over
 *      fix with much tighter targeting and near-zero cost when it fails.
 *   6. Replacing the exact solver's line-centric include/exclude tree with a
 *      coverage-canonical point-pivot search was supposed to cut state
 *      explosion by branching on one uncovered point and encoding the no-heavy
 *      case with point-level residual bans instead of blocked heavy-line IDs.
 *      Three 770 attempts all failed much earlier than baseline: the
 *      active-mask-only version timed out at N=740 after 30.051s with 154475
 *      nodes, the point-barred exact version also timed out at N=740 after
 *      30.054s with 153245 nodes, and the gain-weighted pivot follow-up still
 *      timed out there at 30.184s with 672732 nodes, versus the 8.672198s
 *      baseline at N=740. The simplest explanation is that these point pivots
 *      mostly branch on many low-gain local choices and lose the strong
 *      pruning and ordering effect of the original high-gain line-centric
 *      search, so the tree gets much deeper before the Lagrangian bound can
 *      bite. Future attempts should avoid replacing heavy-line-first branching
 *      with point-centric fan-out unless a much tighter point-branch upper
 *      bound or forced-reduction scheme is available to keep low-gain pivot
 *      choices from dominating the search.
 *   7. Dynamic and static scarcity-based branch tie-breaks were supposed to
 *      speed up the exact search by steering equal-gain choices toward lines
 *      covering points with fewer alternatives. The live productive-degree
 *      version finished the 770 sweep at 80.377s and 79.327s versus the
 *      72.916s and 73.139s baseline, and three cheaper static follow-ups still
 *      lost at 124.369s, 73.842s, and 77.868s with the hard N=758-N=769 tail
 *      never beating baseline. The simplest explanation is that scarce-point
 *      ordering either adds bookkeeping cost or pulls the proof into deeper
 *      local trees without lifting the raw late-tail bound, so future attempts
 *      should avoid branch heuristics that promote local rarity unless they
 *      can be shown to strengthen pruning on the same hard tail cases.
 *   8. Duplicate-state pruning via StateKey/seen_search_state() and the
 *      matching frontier cache was suspected of hurting performance because it
 *      rebuilds and hashes the active mask plus a filtered blocked-line set at
 *      many search states. A full removal lost at 76.024s, the DFS-only
 *      removal lost at 75.986s, the blocked-line gated version lost at
 *      75.677s, and the frontier-only variant split 71.931s / 73.348s against
 *      the 72.849s / 71.992s baseline pair, so none beat baseline
 *      repeatably. The simplest explanation is that the duplicate hits still
 *      save enough repeated search to pay for their bookkeeping, so this is
 *      likely not the bottleneck. Future attempts should avoid weakening state
 *      dedup unless they preserve most of its pruning or replace it with a
 *      demonstrably cheaper equivalent signature on the same hard tail cases.
 *   9. The pressure-based branch tie-break and its ordered_point_lines_
 *      support structure were suspected of hurting performance because they
 *      rebuild per-point helper data for every N and spend extra work inside
 *      choose_branch_line() on a score that is not fully current-path aware.
 *      In practice the full removal regressed badly to 84.797s / 83.635s
 *      against the 73.065s / 68.579s baseline pair, the dead-sort-only follow-up
 *      only managed a near-tie at 70.668s / 70.123s, the raw-incidence proxy
 *      landed at 70.971s, and the thresholded version fell back to 82.689s.
 *      The simplest explanation is that the pressure heuristic still saves
 *      enough late-tail search to outweigh its bookkeeping, while the unused
 *      sort itself is only a small constant-cost leak rather than the true
 *      bottleneck. Future attempts should avoid weakening this branch-ordering
 *      logic unless the replacement keeps the same hard-tail pruning power and
 *      removes substantially more work than a small setup-only cleanup.
 *   10. Branching on minimum-conflict points was supposed to speed up the
 *      exact search by resolving the most constrained active points before the
 *      wider heavy-line tree could branch. In practice the pure MRV version
 *      timed out at N=740 after 30.042s with 102808 nodes versus the
 *      5.980834s / 5.968596s baseline, the gain-first MRV tie-break slowed
 *      N=740/N=742/N=744 to 6.876605s/3.265492s/3.969776s versus the
 *      5.980834s/2.665762s/3.063853s baseline run, and the degree-1-only
 *      fallback still drove N=740 past 17.6s before it was cut. The simplest
 *      explanation is that even selective point-scarcity branching disrupts
 *      the solver's stronger high-gain line order and sends the proof into
 *      deeper local trees before the Lagrangian bound can bite, so future
 *      attempts should avoid MRV-style point steering unless it preserves the
 *      same high-gain branch order and adds almost no extra branch-selection
 *      work.
 *   11. Retiring active points that had no remaining productive heavy-line
 *      coverage was supposed to speed up the DFS by shrinking the active
 *      masks, frontier states, and duplicate-state signatures before the hard
 *      tail search. Three variations all failed: the full recursive version
 *      immediately returned wrong answers (for example N=740 dropped from the
 *      99-line baseline to 94 or 102 lines), the line-cover-adjusted repair
 *      still changed certified results, and the deeper-only version kept the
 *      same wrong early answers while running slower. The simplest explanation
 *      is that these "residual-only" points still participate in the solver's
 *      exact gain/state invariants through line-cover, seeding, and
 *      branch-state interactions in ways this reduction did not preserve.
 *      Future attempts should avoid stripping active points out of the exact
 *      state unless every downstream gain, cover-count, undo, and witness
 *      invariant is rederived around a formally exact reduced representation.
 *   12. Rebuilding lagrangian_upper_bound()'s active-point list from the
 *      current productive-line set was supposed to remove repeated per-point
 *      incidence scans inside the solver's hottest bound routine. The direct
 *      productive-degree version broke exactness immediately (for example
 *      N=740 rose from the 99-line baseline to 102 lines), and two exact
 *      follow-ups that stamped productive-line points globally or only below a
 *      2048-line threshold still lost to the 45.410s / 46.229s baseline pair,
 *      finishing the full 770 sweep in 48.972s and 48.870s with the
 *      N=759-N=769 tail slower. The simplest explanation is that walking the
 *      productive lines and maintaining the extra stamp state costs more than
 *      the original incidence early-exit scan saves on these instances, so
 *      future attempts should avoid replacing the bound's point-discovery step
 *      unless the new path is both provably exact and cheaper on the large-root
 *      and hard-tail cases that dominate runtime.
 *   13. Caching each line's branch-pressure score inside the task state was
 *      supposed to speed up choose_branch_line() by replacing repeated
 *      point-rescans with incremental updates. Two full 770 runs came in at
 *      45.585s / 45.698s versus the 45.147s / 45.702s baseline, the
 *      available-only follow-up immediately regressed N=741 from about 0.058s
 *      to 2.716s with 9047 nodes, and the tail-only gated variant was already
 *      slower through N=745 at 13.019s versus the 11.925s sanity baseline.
 *      The simplest explanation is that maintaining or selectively rebuilding
 *      this extra per-line state costs at least as much as the rescans it
 *      replaces, so future attempts should avoid branch-order caches unless
 *      they can be updated with near-zero copy/undo overhead and proven not to
 *      disturb formerly trivial witness-hit cases.
 *   14. Conflict-graph decomposition on disconnected productive heavy-line
 *       components was supposed to speed up the hard tail by solving
 *       independent subproblems separately instead of letting the DFS branch
 *       across their Cartesian product. Three 770 follow-ups all lost to the
 *       44.470s / 44.238s baseline: the direct nested-solver version regressed
 *       to 64.993s and even broke N=763, the conservative incumbent-seeding
 *       version still took 47.031s, and the tiny exact-component brute-force
 *       version still took 48.124s with the same N=763-N=765 slowdown. The
 *       simplest explanation is that disconnected productive components are
 *       too rare or too small here, so the extra detection and subproblem
 *       bookkeeping costs more than the independence it uncovers. Future
 *       attempts should avoid residual-component decomposition unless a much
 *       cheaper trigger can prove that large independent pockets appear
 *       frequently on the actual hard-tail search states.
 *   15. Adding a remaining-savings pruning bound was supposed to speed up the
 *       exact search by cheaply rejecting branches whose current heavy-line
 *       covers could not recover enough gain. In three variations it either
 *       timed out at N=740 after 30.142s with 120708 nodes, regressed the
 *       full 770 sweep to 47.958s, or split 43.639s / 50.354s against the
 *       45.665s / 44.393s baseline pair, so it did not beat baseline
 *       repeatably. The simplest explanation is that this extra branch-level
 *       savings scan is too loose to prune consistently and still costs enough
 *       on the hard tail to hurt stability. Future attempts should avoid broad
 *       remaining-gain caps unless they are much tighter on the N=763-N=769
 *       cases and cheap enough that a miss adds almost no extra work.
 *   16. Carrying a separate exact productive-point state for
 *       lagrangian_upper_bound() was supposed to cut the per-node bound cost
 *       by updating only the touched points and lines across include,
 *       exclude, and undo. In testing, the first version timed out at N=740
 *       after 70.055s with 134465 nodes, and two follow-ups that kept a
 *       separate lagrangian degree plus either an incremental point list or a
 *       rebuilt point list ran faster overall at 116.475s and 110.511s but
 *       changed certified answers starting at N=742 by dropping from the
 *       100-line baseline to 99. The simplest explanation is that even a
 *       small mismatch between the carried bound-only state and the true
 *       residual state can understate the upper bound, so future attempts
 *       should avoid shadow productive-point state unless every include,
 *       exclude, frontier-copy, and undo path is proven to preserve exactly
 *       the same bound inputs as a full rebuild.
 *   17. Reusing productive-degree point activity and two parent-dual include
 *       caps were supposed to speed up lagrangian_upper_bound() by avoiding
 *       its per-node active-point rediscovery and pruning include children
 *       from carried parent state instead of rebuilding the bound immediately.
 *       In testing, the productive-degree-only version changed certified
 *       answers at N=740 from the 99-line baseline to 102 lines, the direct
 *       include-cap variant also broke N=740, and the cached-dual follow-up
 *       still broke N=742 by rising from 99 lines to 100 lines. The simplest
 *       explanation is that these carried summaries do not exactly match the
 *       full rebuilt bound inputs after seeding, forced reductions, and dual
 *       polishing, so future attempts should avoid pruning or point-activity
 *       reuse from shadow bound state unless it is proven bit-for-bit
 *       equivalent to the rebuilt path.
 *   18. Lowering the forced-line kernel from 4 uniquely owned active points to
 *       selective 3-point cases was supposed to shrink the exact DFS earlier
 *       by forcing more heavy lines before branching. The broad 3-point
 *       version and the exact-3-point-only follow-up both broke exactness as
 *       early as N=702 by returning 92 lines instead of the 91-line baseline,
 *       and the hard-tail-only gated version kept the answers but slowed the
 *       full 780 sweep to 97.478s and 197716 nodes versus the 95.976s /
 *       95.699s and 187618 / 188266-node baseline pair. The simplest
 *       explanation is that 3-point "uniqueness" in this solver is not a safe
 *       forcing criterion once residual pair covers and the exact-gain search
 *       state interact, so future attempts should avoid forcing 3-unique heavy
 *       lines unless that proof is extended to the residual-cover model and
 *       the solver's exact invariants.
 *   19. Rebuilding lagrangian_upper_bound()'s active-point list from current
 *       productive-line metadata was supposed to remove the per-point
 *       incidence scans at the start of the bound. The direct
 *       productive-degree version immediately broke exactness at N=700 by
 *       rising from the 91-line baseline to 92, the exact point-stamp rebuild
 *       finished the 780 sweep at 96.392s with 186658 nodes, and the
 *       productive-mask union follow-up slowed further to 96.937s with 187673
 *       nodes, against the 95.466s / 96.385s and 188978 / 188166-node
 *       baseline pair. The simplest explanation is that the maintained
 *       productive summaries are either not exact enough for pruning or cost
 *       at least as much to rebuild from as the original incidence early-exit
 *       scan, so future attempts should avoid replacing this bound setup
 *       unless the alternative representation is both bit-for-bit exact and
 *       cheaper on the real N=763-N=775 tail states.
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
constexpr int    kStartN               = 700;   // N to start from in the timed sweep; set to 0 to run from 1
constexpr int    kExecutionLimit       = 770;   // hard cap and default run-to N
// Maximum wall-clock seconds allowed per N. Set to 0.0 to disable (no limit). 60=1 min
constexpr double kPerNTimeLimitSeconds = 70;

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
        return static_cast<std::size_t>(LineKeyHash::mix(packed));
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

[[nodiscard]] SlopeKey reduced_slope(int dx, int dy) noexcept {
    const int g = std::gcd(dx, dy);
    return {dy / g, dx / g};
}

[[nodiscard]] bool anchor_is_leftmost_for_slope(
    const std::vector<int>& primes, int anchor, const SlopeKey& slope) noexcept
{
    const int step_x = slope.dx;
    const std::int64_t step_y = slope.dy;
    std::int64_t expected_y =
        static_cast<std::int64_t>(primes[static_cast<std::size_t>(anchor)]) - step_y;
    for (int x = anchor - step_x; x >= 0; x -= step_x, expected_y -= step_y) {
        if (expected_y ==
            static_cast<std::int64_t>(primes[static_cast<std::size_t>(x)]))
            return false;
    }
    return true;
}

[[nodiscard]] std::vector<HeavyLine> enumerate_heavy_lines(const std::vector<int>& primes) {
    const int n = static_cast<int>(primes.size());
    std::vector<HeavyLine> lines;
    lines.reserve(static_cast<std::size_t>(n));

    std::unordered_map<SlopeKey, std::vector<int>, SlopeKeyHash> slope_buckets;
    slope_buckets.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        slope_buckets.clear();
        const int base_prime = primes[static_cast<std::size_t>(i)];
        for (int j = i + 1; j < n; ++j) {
            const SlopeKey slope = reduced_slope(
                j - i,
                primes[static_cast<std::size_t>(j)] - base_prime);
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

struct WitnessCover {
    std::vector<LineKey> lines;
};

[[nodiscard]] bool line_contains_point(const LineKey& line, int x, int y) noexcept {
    return line.a * static_cast<std::int64_t>(x) +
           line.b * static_cast<std::int64_t>(y) +
           line.c == 0;
}

[[nodiscard]] WitnessCover build_witness_cover(
    int n,
    const std::vector<int>& primes,
    const std::vector<HeavyLine>& lines,
    const std::vector<int>& chosen_heavy_line_ids)
{
    WitnessCover witness;
    witness.lines.reserve(
        chosen_heavy_line_ids.size() + static_cast<std::size_t>(ceil_half(n)));

    BitMask1024 covered;
    for (int id : chosen_heavy_line_ids) {
        const auto& heavy = lines[static_cast<std::size_t>(id)];
        witness.lines.push_back(canonical_line(
            heavy.points[0], primes[static_cast<std::size_t>(heavy.points[0])],
            heavy.points[1], primes[static_cast<std::size_t>(heavy.points[1])]));
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
            pending, primes[static_cast<std::size_t>(pending)],
            p,       primes[static_cast<std::size_t>(p)]));
        pending = -1;
    }
    if (pending >= 0)
        witness.lines.push_back({1, 0, -static_cast<std::int64_t>(pending)});

    return witness;
}

[[nodiscard]] bool witness_covers_point(
    const WitnessCover& witness, int x, int y) noexcept
{
    for (const LineKey& line : witness.lines)
        if (line_contains_point(line, x, y))
            return true;
    return false;
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
        apply_forced_line_reduction(root, root_ctx, false);

        const int greedy_gain = greedy_completion(root, root_ctx, true);
        if (greedy_gain > current_best_gain())
            submit_candidate(greedy_gain, root_ctx.greedy_lines);
        greedy_complete_from_seed(root, root_ctx, warm_start_line_ids);

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
        std::vector<int>              productive_degree;
        std::vector<int>              sole_productive_line;
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
        std::size_t covered_point_metadata = 0;
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
        std::vector<int>           forced_unique_count;
        std::vector<int>           forced_touched_lines;
        std::vector<int>           covered_point_degree;
        std::vector<int>           covered_point_sole;
    };

    static constexpr int kNoProductiveLine = -1;
    static constexpr int kManyProductiveLines = -2;

    [[nodiscard]] Task make_root_task(
        int n, const std::vector<int>& active_ids) const
    {
        Task root;
        root.active = BitMask1024::prefix(n);
        root.active_count = n;
        root.line_cover.assign(lines_.size(), 0);
        root.productive_degree.assign(static_cast<std::size_t>(n), 0);
        root.sole_productive_line.assign(
            static_cast<std::size_t>(n), kNoProductiveLine);
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
        rebuild_productive_point_metadata(root);
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

        if (ctx.forced_unique_count.size() != lines_.size())
            ctx.forced_unique_count.resize(lines_.size());
        if (ctx.forced_touched_lines.capacity() < lines_.size())
            ctx.forced_touched_lines.reserve(lines_.size());
        ctx.forced_touched_lines.clear();

        if (ctx.covered_point_degree.capacity() < static_cast<std::size_t>(cur_n_))
            ctx.covered_point_degree.reserve(static_cast<std::size_t>(cur_n_));
        ctx.covered_point_degree.clear();

        if (ctx.covered_point_sole.capacity() < static_cast<std::size_t>(cur_n_))
            ctx.covered_point_sole.reserve(static_cast<std::size_t>(cur_n_));
        ctx.covered_point_sole.clear();

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
            ctx.covered_point_degree.size(),
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

    void rebuild_productive_point_metadata(Task& task) const {
        std::fill(task.productive_degree.begin(), task.productive_degree.end(), 0);
        std::fill(
            task.sole_productive_line.begin(),
            task.sole_productive_line.end(),
            kNoProductiveLine);
        for (int id : task.productive_ids) {
            for (int p : lines_[static_cast<std::size_t>(id)].points) {
                if (p >= cur_n_) break;
                if (!task.active.test(p)) continue;
                int& degree = task.productive_degree[static_cast<std::size_t>(p)];
                int& sole = task.sole_productive_line[static_cast<std::size_t>(p)];
                if (degree == 0) sole = id;
                else             sole = kManyProductiveLines;
                ++degree;
            }
        }
    }

    [[nodiscard]] int find_any_productive_line_for_point(
        const Task& task, int p) const
    {
        for (int id : incidence_[static_cast<std::size_t>(p)]) {
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;
            if (task.productive_pos[static_cast<std::size_t>(id)] >= 0)
                return id;
        }
        return kNoProductiveLine;
    }

    void on_productive_line_removed(Task& task, int id) const {
        for (int p : lines_[static_cast<std::size_t>(id)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;

            int& degree = task.productive_degree[static_cast<std::size_t>(p)];
            int& sole = task.sole_productive_line[static_cast<std::size_t>(p)];
            if (degree <= 0) continue;

            if (degree == 1) {
                degree = 0;
                sole = kNoProductiveLine;
            } else if (degree == 2) {
                degree = 1;
                sole = find_any_productive_line_for_point(task, p);
            } else {
                --degree;
                sole = kManyProductiveLines;
            }
        }
    }

    void on_productive_line_reinserted(Task& task, int id) const {
        for (int p : lines_[static_cast<std::size_t>(id)].points) {
            if (p >= cur_n_) break;
            if (!task.active.test(p)) continue;

            int& degree = task.productive_degree[static_cast<std::size_t>(p)];
            int& sole = task.sole_productive_line[static_cast<std::size_t>(p)];
            if (degree == 0) {
                degree = 1;
                sole = id;
            } else {
                ++degree;
                sole = kManyProductiveLines;
            }
        }
    }

    void remove_productive_impl(
        Task& task, int id, std::vector<int>* removed_log) const
    {
        int& pos = task.productive_pos[static_cast<std::size_t>(id)];
        if (pos < 0) return;
        const int last_id = task.productive_ids.back();
        task.productive_ids[static_cast<std::size_t>(pos)] = last_id;
        task.productive_pos[static_cast<std::size_t>(last_id)] = pos;
        task.productive_ids.pop_back();
        pos = -1;
        if (removed_log != nullptr)
            removed_log->push_back(id);
        on_productive_line_removed(task, id);
    }

    // O(1) swap-to-back removal from productive_ids; records in the worker stack for reversal.
    void remove_productive(Task& task, WorkerContext& ctx, int id) const {
        remove_productive_impl(task, id, &ctx.removed_from_productive);
    }

    // Called during undo to re-insert a line that was removed from the productive set.
    void reinsert_productive(Task& task, int id) const {
        if (task.productive_pos[static_cast<std::size_t>(id)] >= 0) return;
        task.productive_pos[static_cast<std::size_t>(id)] =
            static_cast<int>(task.productive_ids.size());
        task.productive_ids.push_back(id);
        on_productive_line_reinserted(task, id);
    }

    void disable_line_without_undo(Task& task, int id) const {
        unsigned char& available = task.available[static_cast<std::size_t>(id)];
        if (!available) return;
        available = 0;
        task.blocked_line_ids.push_back(id);
        remove_productive_impl(task, id, nullptr);
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
        ctx.covered_point_degree.push_back(
            task.productive_degree[static_cast<std::size_t>(p)]);
        ctx.covered_point_sole.push_back(
            task.sole_productive_line[static_cast<std::size_t>(p)]);
        task.productive_degree[static_cast<std::size_t>(p)] = 0;
        task.sole_productive_line[static_cast<std::size_t>(p)] = kNoProductiveLine;
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
        for (std::size_t i = mark.removed_from_productive;
             i < ctx.removed_from_productive.size(); ++i)
            reinsert_productive(task, ctx.removed_from_productive[i]);
        for (std::size_t i = mark.covered_point_metadata;
             i < ctx.covered_point_degree.size(); ++i) {
            const int p = ctx.newly_covered_points[i];
            task.productive_degree[static_cast<std::size_t>(p)] =
                ctx.covered_point_degree[i];
            task.sole_productive_line[static_cast<std::size_t>(p)] =
                ctx.covered_point_sole[i];
        }
        task.blocked_line_ids.resize(mark.blocked_lines);

        ctx.removed_from_productive.resize(mark.removed_from_productive);
        ctx.decremented_lines.resize(mark.decremented_lines);
        ctx.newly_covered_points.resize(mark.newly_covered_points);
        ctx.covered_point_degree.resize(mark.covered_point_metadata);
        ctx.covered_point_sole.resize(mark.covered_point_metadata);
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
                if (old_cover == 3)   // cover 3->2: drop from productive set
                    remove_productive_impl(task, lid, nullptr);
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

    void greedy_complete_from_seed(
        const Task& base,
        WorkerContext& ctx,
        const std::vector<int>& line_ids)
    {
        if (line_ids.empty()) return;

        Task trial = base;
        ensure_task_search_buffers(trial);
        for (int id : line_ids) {
            if (id < 0 || id >= static_cast<int>(lines_.size())) continue;
            if (lines_[static_cast<std::size_t>(id)].activate_at > cur_n_) continue;

            const int delta = positive_gain(trial, id);
            if (delta <= 0) continue;

            select_line_without_undo(trial, id);
            trial.current_gain += delta;
            trial.current_choice.push_back(id);
        }
        maybe_submit_greedy(trial, ctx);
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

    void apply_forced_line_reduction(
        Task& task, WorkerContext& ctx, bool record_undo) const
    {
        auto& unique_count = ctx.forced_unique_count;
        auto& touched = ctx.forced_touched_lines;

        while (true) {
            touched.clear();

            task.active.for_each_set_bit(words_, [&](int p) {
                const int unique_line =
                    task.sole_productive_line[static_cast<std::size_t>(p)];
                if (unique_line >= 0) {
                    int& count = unique_count[static_cast<std::size_t>(unique_line)];
                    if (count++ == 0)
                        touched.push_back(unique_line);
                }
            });

            int best_line = -1;
            int best_unique = 2;
            int best_gain = -1;
            for (int id : touched) {
                const int unique_points = unique_count[static_cast<std::size_t>(id)];
                unique_count[static_cast<std::size_t>(id)] = 0;
                if (unique_points < 4) continue;

                const int gain = positive_gain(task, id);
                if (gain <= 0) continue;

                const bool better =
                    (unique_points > best_unique) ||
                    (unique_points == best_unique && gain > best_gain) ||
                    (unique_points == best_unique && gain == best_gain &&
                     (best_line < 0 || id < best_line));
                if (better) {
                    best_line = id;
                    best_unique = unique_points;
                    best_gain = gain;
                }
            }

            if (best_line < 0) break;

            // Safe kernel rule: if a productive heavy line owns at least three
            // currently active points that no other productive heavy line can
            // cover, there is always an optimal completion that includes it.
            if (record_undo) select_line(task, ctx, best_line);
            else             select_line_without_undo(task, best_line);
            task.current_gain += best_gain;
            task.current_choice.push_back(best_line);
        }
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
                        std::vector<long double>* subgradient) -> long double {
            long double value = 0.0L;

            for (int p : active_points) {
                value += yy[static_cast<std::size_t>(p)];
                if (subgradient != nullptr)
                    (*subgradient)[static_cast<std::size_t>(p)] = 1.0L;
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
                    if (subgradient != nullptr) {
                        for (int p : lines_[static_cast<std::size_t>(id)].points) {
                            if (p >= cur_n_) break;
                            if (task.active.test(p))
                                (*subgradient)[static_cast<std::size_t>(p)] -= 1.0L;
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

        constexpr int kCoordinatePolishMaxSweeps = 6;
        if (kCoordinatePolishMaxSweeps > 0) {
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
            for (int sweep = 0; sweep < kCoordinatePolishMaxSweeps; ++sweep) {
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
                if (gain_upper_bound_cannot_beat_best_cost(
                        static_cast<long double>(task.current_gain) + polished_best))
                    return static_cast<long double>(task.current_gain) + polished_best;
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
        const UndoMark node_undo = mark_undo(task, ctx);
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
            rollback();
            return;
        }
        const long double ub = lagrangian_upper_bound(task, ctx);
        capture_root_cost_lb(ub);
        if (root_interval_closed()) {
            rollback();
            return;
        }
        if (gain_upper_bound_cannot_beat_best_cost(ub)) {
            rollback();
            return;
        }

        maybe_submit_greedy(task, ctx);
        if (root_interval_closed()) {
            rollback();
            return;
        }

        const int branch_line = choose_branch_line(task);
        if (branch_line < 0) {
            submit_candidate(task.current_gain, task.current_choice);
            rollback();
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

        if (gain_upper_bound_cannot_beat_best_cost(exclude_ub)) {
            rollback();
            return;
        }

        {
            const UndoMark undo = mark_undo(task, ctx);
            disable_line(task, ctx, branch_line);
            dfs(task, ctx);
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
        frontier.reserve(std::max<std::size_t>(
            4, static_cast<std::size_t>(worker_count_) * 8U));
        frontier.push_back(std::move(root));

        std::vector<Task> children;
        const unsigned frontier_multiplier =
            current_best_cost() >= 93 ? 8U : 4U;
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
            apply_forced_line_reduction(task, ctx, false);
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
                apply_forced_line_reduction(child, ctx, false);
                if (remember_frontier_state(child, ctx, frontier_best))
                    continue;
                frontier.push_back(std::move(child));
            }
        }

        std::sort(frontier.begin(), frontier.end(), [](const Task& a, const Task& b) {
            if (a.active_count != b.active_count) return a.active_count > b.active_count;
            if (a.productive_ids.size() != b.productive_ids.size())
                return a.productive_ids.size() > b.productive_ids.size();
            return a.current_gain < b.current_gain;
        });

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
struct CoverInequalityResult {
    int                        lb = 0;
    std::vector<unsigned char> in_support;
};

[[nodiscard]] CoverInequalityResult cover_inequality_lb(
    int n, const std::vector<int>& active_ids, const std::vector<HeavyLine>& lines)
{
    CoverInequalityResult result;
    if (n == 0) return result;

    result.in_support.assign(static_cast<std::size_t>(n), 1);
    if (active_ids.empty()) {
        result.lb = (n + 1) / 2;
        return result;
    }

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
            if (result.in_support[static_cast<std::size_t>(p)] &&
                violation_count[static_cast<std::size_t>(p)] > worst_v) {
                worst_v = violation_count[static_cast<std::size_t>(p)];
                worst_p = p;
            }
        }
        if (worst_v == 0) break;

        result.in_support[static_cast<std::size_t>(worst_p)] = 0;
        --S_size;

        for (int i : local_inc[static_cast<std::size_t>(worst_p)]) {
            const bool was_violated = (line_size[i] >= 3);
            --line_size[i];
            if (was_violated && line_size[i] < 3) {
                // Line i just fell below the threshold; un-charge its S members.
                for (int q : lines[static_cast<std::size_t>(active_ids[i])].points) {
                    if (q >= n) break;
                    if (result.in_support[static_cast<std::size_t>(q)])
                        --violation_count[static_cast<std::size_t>(q)];
                }
            }
        }
    }

    result.lb = (S_size + 1) / 2;  // ceil(|S| / 2)
    return result;
}

// Root-only LP-relaxation lower bound for
//   min heavy_lines + 0.5 * residual_points.
//
// Its dual is
//   max sum_p y_p
//   s.t. sum_{p in line} y_p <= 1 for every heavy line,
//        0 <= y_p <= 0.5.
//
// Any feasible dual gives a valid lower bound, and ceil(dual_value) remains
// valid because the true optimum is integral.
[[nodiscard]] int fractional_cover_lp_lb(
    int n,
    const std::vector<int>& active_ids,
    const std::vector<HeavyLine>& lines,
    const std::vector<std::vector<int>>& incidence,
    const std::vector<unsigned char>& cover_support)
{
    if (n == 0) return 0;

    constexpr long double kPointCap = 0.5L;
    constexpr long double kLineCap = 1.0L;
    constexpr long double kEps = 1e-18L;

    std::vector<unsigned char> active_line(lines.size(), 0);
    for (int id : active_ids)
        active_line[static_cast<std::size_t>(id)] = 1;

    std::vector<int> degree(static_cast<std::size_t>(n), 0);
    for (int p = 0; p < n; ++p)
        for (int id : incidence[static_cast<std::size_t>(p)])
            if (active_line[static_cast<std::size_t>(id)])
                ++degree[static_cast<std::size_t>(p)];

    std::vector<int> point_order(static_cast<std::size_t>(n));
    std::iota(point_order.begin(), point_order.end(), 0);
    std::stable_sort(point_order.begin(), point_order.end(), [&](int lhs, int rhs) {
        const int dl = degree[static_cast<std::size_t>(lhs)];
        const int dr = degree[static_cast<std::size_t>(rhs)];
        if (dl != dr) return dl < dr;
        return lhs < rhs;
    });

    auto improve_feasible_dual = [&](std::vector<long double> y) -> long double {
        std::vector<long double> line_sum(lines.size(), 0.0L);
        long double total = 0.0L;
        for (int p = 0; p < n; ++p) {
            const long double yp = y[static_cast<std::size_t>(p)];
            if (yp <= kEps) continue;
            total += yp;
            for (int id : incidence[static_cast<std::size_t>(p)]) {
                if (!active_line[static_cast<std::size_t>(id)]) continue;
                line_sum[static_cast<std::size_t>(id)] += yp;
            }
        }

        auto add_max_to_point = [&](int p) -> long double {
            long double cap = kPointCap - y[static_cast<std::size_t>(p)];
            if (cap <= kEps) return 0.0L;
            for (int id : incidence[static_cast<std::size_t>(p)]) {
                if (!active_line[static_cast<std::size_t>(id)]) continue;
                cap = std::min(
                    cap, kLineCap - line_sum[static_cast<std::size_t>(id)]);
                if (cap <= kEps) return 0.0L;
            }
            y[static_cast<std::size_t>(p)] += cap;
            total += cap;
            for (int id : incidence[static_cast<std::size_t>(p)]) {
                if (!active_line[static_cast<std::size_t>(id)]) continue;
                line_sum[static_cast<std::size_t>(id)] += cap;
            }
            return cap;
        };

        for (int sweep = 0; sweep < 3; ++sweep) {
            bool changed = false;
            for (int p : point_order)
                if (add_max_to_point(p) > kEps)
                    changed = true;
            if (!changed) break;
        }
        return total;
    };

    std::vector<long double> seeded_from_cover(static_cast<std::size_t>(n), 0.0L);
    for (int p = 0; p < n; ++p) {
        if (p < static_cast<int>(cover_support.size()) &&
            cover_support[static_cast<std::size_t>(p)]) {
            seeded_from_cover[static_cast<std::size_t>(p)] = kPointCap;
        }
    }

    const long double cover_value = improve_feasible_dual(seeded_from_cover);
    const long double sparse_value =
        improve_feasible_dual(std::vector<long double>(static_cast<std::size_t>(n), 0.0L));
    const long double best_value = std::max(cover_value, sparse_value);
    return static_cast<int>(std::ceil(best_value - 1e-12L));
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
    const CoverInequalityResult cover = cover_inequality_lb(n, active_ids, lines);
    const int fractional_lb = fractional_cover_lp_lb(
        n, active_ids, lines, incidence, cover.in_support);
    const int effective_lb_floor = std::max({root_lb_floor, cover.lb, fractional_lb});
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
        : program_start_(program_start)
        , thread_([this] { run(); })
    {}

    ~LiveDisplay() {
        active_.store(false, std::memory_order_relaxed);
        if (thread_.joinable()) thread_.join();
        clear_line();
    }

    void begin_n(int n) noexcept {
        n_start_ns_.store(
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - program_start_).count(),
            std::memory_order_relaxed);
        current_n_.store(n, std::memory_order_relaxed);
    }

    void clear() const noexcept {
        current_n_.store(0, std::memory_order_relaxed);
        n_start_ns_.store(0, std::memory_order_relaxed);
        clear_line();
    }

private:
    static constexpr auto   kRefreshPeriod = std::chrono::milliseconds(100);
    static constexpr double kActivationDelaySeconds = 0.2;

    void clear_line() const {
        static constexpr char kBlankLine[] =
            "\r                                                                              \r";
        std::cerr << kBlankLine << std::flush;
    }

    void run() {
        while (active_.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(kRefreshPeriod);
            if (!active_.load(std::memory_order_relaxed)) break;

            const int n = current_n_.load(std::memory_order_relaxed);
            if (n == 0) continue;

            const std::int64_t start_ns = n_start_ns_.load(std::memory_order_relaxed);
            if (start_ns == 0) continue;

            const std::int64_t now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - program_start_).count();
            const double n_elapsed = static_cast<double>(now_ns - start_ns) * 1e-9;
            if (n_elapsed < kActivationDelaySeconds) continue;

            const double total = static_cast<double>(now_ns) * 1e-9;
            char buffer[96];
            const int written = std::snprintf(
                buffer, sizeof(buffer),
                "\r[running] N=%d  N_time=%.1fs  total=%.1fs  ",
                n, n_elapsed, total);
            if (written > 0)
                std::cerr.write(
                    buffer,
                    static_cast<std::streamsize>(std::min<int>(
                        written, static_cast<int>(sizeof(buffer) - 1))))
                    << std::flush;
        }
    }

    Clock::time_point         program_start_;
    std::atomic<bool>         active_{true};
    mutable std::atomic<int>  current_n_{0};
    mutable std::atomic<std::int64_t> n_start_ns_{0};
    std::thread               thread_;
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

    const auto program_start = std::chrono::steady_clock::now();
    LiveDisplay display(program_start);

    int  prev    = (kStartN > 1) ? kStartN - 1 : 0;
    long long total_nodes_processed = 0;
    std::vector<int> warm_start_solution;
    WitnessCover witness_cover;
    bool witness_ready = false;
    for (int n_init = 1; n_init < kStartN; ++n_init)
        for (int id : activation[static_cast<std::size_t>(n_init)])
            active_ids.push_back(id);
    for (int n = kStartN; n <= requested_n; ++n) {
        for (int id : activation[static_cast<std::size_t>(n)])
            active_ids.push_back(id);

        display.begin_n(n);
        const auto t0 = std::chrono::steady_clock::now();

        const GreedyCoverSolution greedy = greedy_upper_bound_solution(n, active_ids, lines);
        const int ub_raw = greedy.cost;
        int incumbent = std::min(n, prev + 1);
        incumbent     = std::min(incumbent, greedy.cost);
        const int lb0_floor = warm_start_solution.empty() ? 0 : prev;
        const std::size_t active_count = active_ids.size();
        const bool witness_hit =
            witness_ready &&
            witness_covers_point(
                witness_cover,
                n - 1,
                primes[static_cast<std::size_t>(n - 1)]);
        if (witness_hit)
            incumbent = std::min(incumbent, prev);
        const int ub0_val = incumbent;

        bool timed_out = false;
        int out_lb = 0;
        int out_lb_raw = 0;
        long long out_nodes = 0;
        int answer = 0;
        if (witness_hit) {
            answer = prev;
            out_lb = prev;
            out_lb_raw = prev;
        } else {
            answer = solve_exact(
                n, active_ids, lines, incidence, lb0_floor, greedy.line_ids,
                warm_start_solution, timed_out, out_lb, out_lb_raw, out_nodes);
            if (!timed_out) {
                witness_cover =
                    build_witness_cover(n, primes, lines, warm_start_solution);
                witness_ready = true;
            }
        }
        const int gap = ub0_val - out_lb;
        const int gap_raw = ub_raw - out_lb_raw;

        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t0).count();

        display.clear();
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
    return 0;
}