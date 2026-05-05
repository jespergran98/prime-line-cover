#!/usr/bin/env bash
# =============================================================================
#  PRIMECOVER FRONTIER DIAGNOSIS TOOL
# =============================================================================
#
#  WHAT THIS TOOL DOES
#  ───────────────────
#  This script compiles and runs the primecover solver, then displays a live
#  dashboard showing its progress and memory usage. Once the solver reaches
#  N=808, it switches to a safety table that tells you which frontier
#  multipliers (4U → 131072U) are safe to use on your hardware.
#
#  ⚠️  This tool is ONLY necessary if you plan to run the solver for a "serious"
#      attempt (N > ~850, where the frontier multiplier becomes critical).
#      For quick test runs (N ≤ 850) or if you never change the multiplier,
#      you can ignore this diagnostic tool entirely – the default settings
#      are safe for small N and will not consume excessive RAM.
#
#  HOW TO USE IT
#  ─────────────
#  1. ⚠️  IMPORTANT: Edit the SOURCE_PATH below on line 101 – change the
#     example path to the actual location of your primecover1024.cpp file.
#
#  2. Copy paste this entire script directly into your terminal and hit Enter.
#     (Or save it to a file and run it with bash /full/path/to/script.sh)
#
#  3. Watch the live dashboard. When it switches to the safety table, read
#     the "Status" column to find the highest multiplier marked ✓ SAFE.
#     Example: if "8192U" shows ✓ SAFE but "16384U" shows ✗ UNSAFE,
#     then 8192U is your safe maximum.
#
#  4. Set that multiplier in primecover1024.cpp:
#     - Find the "frontier_multiplier" ladder.
#     - Uncomment the line with your safe multiplier (remove // at the start).
#     - Add a colon (:) right before the line that previously started with "current_best_cost() >= 128 ?"
#     - Save, recompile, and run normally.
#     Example (changing from 8192U to 16384U):
#         Before:
#         //: current_best_cost() >= 133 ? 16384U
#             current_best_cost() >= 128 ? 8192U
#         After:
#             current_best_cost() >= 133 ? 16384U
#           : current_best_cost() >= 128 ? 8192U
#
#  🔒 Cleanup: The script automatically kills the solver on Ctrl+C or terminal close.
#      If a stale process remains, run: pkill -f "./primecover"
#
#  WHAT THE SAFETY TABLE MEANS
#  ────────────────────────────
#  The solver's frontier size = (number of CPU threads) × multiplier.
#  A larger multiplier = more parallel tasks = more RAM consumed.
#  The table predicts RAM usage for each multiplier and compares it against
#  your machine's available RAM:
#
#    ✓  SAFE     — predicted usage is well within your available RAM (< 75%)
#    ⚠  MARGINAL — predicted usage is close to your RAM limit (75–100%)
#    ✗  UNSAFE   — predicted usage exceeds your RAM; the solver will likely OOM
#
#  HOW THE PREDICTION WORKS
#  ─────────────────────────
#  Per-task memory is not constant — it shrinks as the frontier grows, because
#  each task represents a shallower subtree. A naive (RSS / frontier) estimate
#  overstates the true cost at large frontier sizes.
#
#  This tool fits a hyperbolic-decay model to live samples:
#
#    per_task(F)  =  a  +  b / F
#                   ─    ──────
#            asymptotic   overhead that vanishes as F → ∞
#            per-task KB
#
#  Rearranged into a linear form for OLS regression:
#
#    Y  =  RSS / F          (observed per-task memory)
#    X  =  1   / F          (inverse frontier size)
#    Y  =  a  +  b · X      ← fitted by OLS
#
#  The intercept  a  is the asymptotic per-task cost (KB/task).
#  Prediction:    RSS_KB = a · F + b
#
#  This extrapolates accurately even from small-frontier samples, because the
#  curvature of per_task vs F is visible and the intercept is well-identified.
#  Samples are only collected from N=808 onward, where memory behaviour is
#  representative of the full run.
#
#  TIP: Letting the solver continue past N=808 to N=851 (or until the frontier
#  reaches ≥ 8 192 tasks) gives the regression more data points and produces a
#  more accurate safety table.
#
# =============================================================================


# =============================================================================
#  CONFIGURATION  –  set your path here before running
# =============================================================================

# Path to your primecover1024.cpp source file.
# WSL users: Windows drives are mounted under /mnt/, so:
#   C:\Users\yourname\project\primecover1024.cpp
#   becomes:  /mnt/c/Users/yourname/project/primecover1024.cpp
SOURCE_PATH="/mnt/c/path/to/primecover1024.cpp"

# =============================================================================


# ── ANSI colours ──────────────────────────────────────────────────────────────
R='\033[0m'
BOLD='\033[1m'
DIM='\033[2m'
BLUE='\033[38;5;75m'
CYAN='\033[38;5;81m'
GREEN='\033[38;5;83m'
YELLOW='\033[38;5;220m'
RED='\033[38;5;203m'
PURPLE='\033[38;5;141m'
GRAY='\033[38;5;245m'
WHITE='\033[38;5;253m'
ORANGE='\033[38;5;215m'

# ── Box-drawing helpers ────────────────────────────────────────────────────────
W=78   # total inner width (excluding the two border chars)
hline() { printf "${BLUE}${BOLD}┌"; printf '─%.0s' $(seq 1 $W); printf "┐${R}\n"; }
mline() { printf "${BLUE}${BOLD}├"; printf '─%.0s' $(seq 1 $W); printf "┤${R}\n"; }
bline() { printf "${BLUE}${BOLD}└"; printf '─%.0s' $(seq 1 $W); printf "┘${R}\n"; }
row()   { printf "${BLUE}${BOLD}│${R}  %-$((W-2))s${BLUE}${BOLD}│${R}\n" "$*"; }
rowf()  {
  # Render the caller's format+args into a variable, then measure the
  # *visual* width (ANSI escape sequences are zero-width) so we can pad
  # exactly to column W-2 before printing the closing border.
  # wc -m counts characters (not bytes) so multibyte glyphs (─ █ ░) count correctly.
  local content
  content=$(printf "$@")
  local visual
  visual=$(printf '%s' "$content" | sed 's/\x1b\[[0-9;]*m//g')
  local vlen
  vlen=$(printf '%s' "$visual" | LC_ALL=C.UTF-8 wc -m | tr -d ' \t')
  local pad=$(( W - 2 - vlen )); [[ $pad -lt 0 ]] && pad=0
  printf "${BLUE}${BOLD}│${R}  %s%*s${BLUE}${BOLD}│${R}\n" "$content" "$pad" ""
}
empty() { row ""; }

# ── Header (printed once) ─────────────────────────────────────────────────────
print_header() {
  clear
  printf "\n"
  # Top border
  printf "${BLUE}${BOLD}╔"; printf '═%.0s' $(seq 1 $W); printf "╗${R}\n"
  # Header line
  local header_str="  PRIMECOVER LIVE DIAGNOSTICS  ·  PID $SOLVER_PID"
  # Compute visible characters (strip ANSI codes first)
  local plain_str=$(printf '%s' "$header_str" | sed 's/\x1b\[[0-9;]*m//g')
  local visual_len=$(LC_ALL=C.UTF-8 printf '%s' "$plain_str" | wc -m | tr -d ' ')
  local pad=$(( W - visual_len ))
  printf "${BLUE}${BOLD}║${R}%s%*s${BLUE}${BOLD}║${R}\n" "$header_str" $pad ""
  # Bottom border
  printf "${BLUE}${BOLD}╚"; printf '═%.0s' $(seq 1 $W); printf "╝${R}\n"
  printf "\n"
}

# ── Progress bar (Unicode blocks, 50-wide) ────────────────────────────────────
make_bar() {
  local pct=$1 width=50
  local filled=$(( pct * width / 100 ))
  local empty_w=$(( width - filled ))
  printf "${GREEN}${BOLD}"
  printf '█%.0s' $(seq 1 $filled 2>/dev/null)
  printf "${DIM}${GRAY}"
  printf '░%.0s' $(seq 1 $empty_w 2>/dev/null)
  printf "${R}"
}

# ── Safety colour ─────────────────────────────────────────────────────────────
safety_color() {
  local est=$1 ram=$2
  if (( $(echo "$est < $ram * 0.75" | bc -l) )); then
    printf "${GREEN}${BOLD}"
  elif (( $(echo "$est < $ram" | bc -l) )); then
    printf "${YELLOW}${BOLD}"
  else
    printf "${RED}${BOLD}"
  fi
}
safety_icon() {
  local est=$1 ram=$2
  if (( $(echo "$est < $ram * 0.75" | bc -l) )); then echo "✓  SAFE    "
  elif (( $(echo "$est < $ram" | bc -l) )); then     echo "⚠  MARGINAL"
  else                                                echo "✗  UNSAFE  "
  fi
}

# ── Human-readable numbers ────────────────────────────────────────────────────
human_kb() {   # KB → "1 234 MB" style
  local kb=$1
  local mb=$(( kb / 1024 ))
  printf "%'d MB" $mb
}
comma() { printf "%'d" "$1"; }

# =============================================================================
#  SAMPLE ACCUMULATION
#  Each call to read_stats appends one (frontier, rss_kb) pair to SAMPLE_FILE.
#  The regression worker deduplicates by frontier and fits RSS = α + β × F.
# =============================================================================
SAMPLE_FILE="/tmp/pc_diag_samples_$$.tsv"

append_sample() {
  # Only record samples from N=808 onward — memory behaviour before that
  # point is structurally different and would bias the regression.
  [[ "${LAST_N:-0}" -ge 808 ]] || return
  [[ "${FRONTIER:-0}" -gt 0 && "${RSS_KB:-0}" -gt 0 ]] || return
  printf '%d\t%d\n' "$FRONTIER" "$RSS_KB" >> "$SAMPLE_FILE"
}

# =============================================================================
#  REGRESSION ENGINE  (hyperbolic-decay model)
#  Reads SAMPLE_FILE, deduplicates (keeping max RSS per frontier bucket so we
#  capture each size's peak pressure), then fits the transformed system:
#
#    Y = RSS / F   (per-task memory, KB/task)
#    X = 1   / F   (inverse frontier)
#    Y = a + b · X   by OLS
#
#  The intercept  a  is the asymptotic per-task cost, independent of how
#  inflated per-task measurements are at small frontiers.
#  Prediction uses  RSS_KB = a · F + b  (same formula shape, correct values).
#
#  Outputs globals:
#    ASYM_READY        – 1 if regression succeeded, 0 otherwise
#    ASYM_BETA_KB      – a in KB/task  (asymptotic per-task estimate)
#    ASYM_ALPHA_MB     – b in MB       (hyperbolic offset, vanishes at large F)
#    ASYM_UNIQUE       – number of distinct frontier sizes used in the fit
#    ASYM_PER_TASK     – best available per-task estimate (a if ready, else PER_TASK)
# =============================================================================
compute_regression() {
  ASYM_READY=0
  ASYM_BETA_KB="0.0"
  ASYM_ALPHA_MB="0"
  ASYM_UNIQUE=0
  ASYM_PER_TASK="$PER_TASK"

  [[ -f "$SAMPLE_FILE" ]] || return

  # awk does the heavy lifting:
  #   1. Deduplicate: for each frontier F, keep the maximum RSS_KB observed.
  #      (Peak RSS is the conservative, correct value to plan against.)
  #   2. Check we have ≥2 distinct frontier buckets (otherwise the system is
  #      under-determined and we cannot separate α from β).
  #   3. Compute OLS:
  #        Σ  →  sf, sr, sf2, sfr  (sums over unique-frontier points)
  #        β  =  (n·Σfr − Σf·Σr) / (n·Σf² − (Σf)²)
  #        α  =  (Σr − β·Σf) / n
  #   4. Sanity-check β > 0 (memory cannot decrease with more tasks).
  local result
  result=$(awk '
    BEGIN { FS="\t" }
    {
      f = $1 + 0; r = $2 + 0
      if (f > 0 && r > 0 && r > max_r[f]) max_r[f] = r
    }
    END {
      # Transform to hyperbolic variables:
      #   Y = per-task memory  = RSS / F
      #   X = inverse frontier = 1   / F
      # Fit Y = a + b*X by OLS.
      # Intercept a = asymptotic per-task KB/task (what we want for predictions).
      # Slope     b = hyperbolic offset  (positive; shrinks as F grows).
      n = 0; sx = 0; sy = 0; sx2 = 0; sxy = 0
      for (f in max_r) {
        r = max_r[f]
        x = 1.0 / f          # X = 1/F
        y = r / (f + 0.0)    # Y = RSS/F  (per-task, KB/task)
        n++; sx += x; sy += y; sx2 += x*x; sxy += x*y
      }
      if (n < 2) { print "INSUFFICIENT " n; exit }
      denom = n * sx2 - sx * sx
      if (denom == 0) { print "SINGULAR"; exit }
      b_slope     = (n * sxy - sx * sy) / denom   # slope  b  (KB)
      a_intercept = (sy - b_slope * sx) / n       # intercept a (KB/task, asymptotic)
      if (a_intercept <= 0) { print "DEGENERATE a=" a_intercept; exit }
      # b_slope can legitimately be negative in noisy data; clamp only for display
      b_mb_display = (b_slope > 0 ? b_slope : 0) / 1024
      printf "OK n=%d beta=%.6f alpha_mb=%.4f b_kb=%.2f\n", \
             n, a_intercept, b_mb_display, b_slope
    }
  ' "$SAMPLE_FILE")

  case "$result" in
    OK*)
      ASYM_UNIQUE=$(  echo "$result" | grep -o 'n=[0-9]*'        | cut -d= -f2)
      ASYM_BETA_KB=$( echo "$result" | grep -o 'beta=[0-9.]*'    | cut -d= -f2)
      ASYM_ALPHA_MB=$(echo "$result" | grep -o 'alpha_mb=[0-9.]*'| cut -d= -f2)
      ASYM_B_KB=$(    echo "$result" | grep -oE 'b_kb=-?[0-9.]+' | cut -d= -f2)
      # With the hyperbolic model, ASYM_BETA_KB (the asymptotic intercept a)
      # is naturally ≤ PER_TASK (since PER_TASK = RSS/F includes the α/F
      # inflation at small frontiers).  No additional sanity gate is needed —
      # the old "beta > 2×PER_TASK" guard was what caused the wild swings when
      # line-count jumps briefly pushed PER_TASK high.
      ASYM_READY=1
      ASYM_PER_TASK="$ASYM_BETA_KB"
      ;;
    *)
      # INSUFFICIENT / SINGULAR / DEGENERATE — leave ASYM_READY=0
      ASYM_UNIQUE=$(echo "$result" | grep -o '[0-9]*$' || echo 0)
      ;;
  esac
}

# =============================================================================
#  RAM PREDICTION
#  Given a task count and the regression result, estimate total RSS in GB.
#  When regression is not ready, falls back to naive scaling.
# =============================================================================
predict_rss_gb() {
  local tasks=$1
  if [[ "$ASYM_READY" -eq 1 ]]; then
    # Hyperbolic model: RSS_KB = a * F + b
    #   a = ASYM_BETA_KB   (asymptotic per-task, KB/task)
    #   b = ASYM_B_KB      (hyperbolic offset, KB)
    #
    # SAFETY CLAMP: b physically represents overhead that vanishes as F->oo and
    # must be >= 0.  Regression noise on sparse data can produce a small negative
    # b_slope, which would cause predict_rss_gb to UNDERESTIMATE RAM usage and
    # incorrectly classify an unsafe multiplier as SAFE or MARGINAL.
    # The display already clamps b to 0 (via alpha_mb).  We must do the same
    # here so the prediction is always at least a*F alone.
    local b_kb_raw="${ASYM_B_KB:-0}"
    local b_kb
    b_kb=$(echo "scale=6; if ($b_kb_raw < 0) 0 else $b_kb_raw" | bc)
    echo "scale=2; ($ASYM_BETA_KB * $tasks + $b_kb) / 1048576" | bc
  else
    # Naive: current per-task × tasks
    echo "scale=2; $tasks * ${PER_TASK:-0} / 1048576" | bc
  fi
}

# ── Read live stats from log ──────────────────────────────────────────────────
read_stats() {
  LAST_N=$(strings primecover.log 2>/dev/null \
    | grep -a 'N=' | tail -1 | sed -n 's/.*N=\([0-9]\+\).*/\1/p')
  LINES=$(strings primecover.log 2>/dev/null \
    | grep -a 'lines=' | tail -1 | sed -n 's/.*lines=\([0-9]\+\).*/\1/p')
  FRONTIER=$(strings primecover.log 2>/dev/null \
    | grep -aE 'frontier=[1-9][0-9]*' | tail -1 \
    | sed -n 's/.*frontier=\([0-9]\+\).*/\1/p')
  [[ -z "$FRONTIER" ]] && FRONTIER=0

  RSS_KB=$(ps -o rss= -p "$SOLVER_PID" 2>/dev/null | tr -d ' ')
  [[ -z "$RSS_KB" ]] && RSS_KB=0

  if [[ "$FRONTIER" -gt 0 && "$RSS_KB" -gt 0 ]]; then
    PER_TASK=$(echo "scale=1; $RSS_KB / $FRONTIER" | bc)
  else
    PER_TASK="0.0"
  fi

  # ── Accumulate regression sample ──────────────────────────────────────────
  append_sample

  # ── Update regression model ───────────────────────────────────────────────
  compute_regression

  VCPUS=$(nproc)
  RAM_GB=$(free -g | awk '/^Mem:/{print $2}')

  if [[ -n "$LAST_N" && "$LAST_N" -ge 0 ]]; then
    # Non-linear mapping: N → fraction of total wall-clock time.
    # Calibrated from observed run: N=500 ≈ 2.5 s (4%), N=780 ≈ 15 s (23%), N=808 ≈ 65 s (100%).
    # Piecewise-linear in time-fraction space so the bar crawls early and
    # accelerates as the hard N values approach – matching solver reality.
    local n=$LAST_N
    if   [[ $n -ge 808 ]]; then PCT=100
    elif [[ $n -ge 780 ]]; then PCT=$(( 23 + (n - 780) * 77 / 28 ))
    elif [[ $n -ge 500 ]]; then PCT=$(( 4  + (n - 500) * 19 / 280 ))
    else                        PCT=$(( n  * 4 / 500 ))
    fi
    [[ $PCT -gt 100 ]] && PCT=100
  else
    PCT=0
  fi
}

# ── Regression status label (one-liner for display boxes) ─────────────────────
regression_status_line() {
  if [[ "$ASYM_READY" -eq 1 ]]; then
    # Format: "a = XX.XX KB/task  (b≈XXXX MB, N pts)"  ≤ 63 visible chars
    local a_short; a_short=$(printf '%.2f' "$ASYM_BETA_KB")
    local b_short; b_short=$(printf '%.0f' "$ASYM_ALPHA_MB")
    printf "${GREEN}${BOLD}a = %s KB/task${R}  ${GRAY}(b ≈ %s MB, %d pts)${R}" \
      "$a_short" "$b_short" "$ASYM_UNIQUE"
  else
    local needed=$(( 2 - ${ASYM_UNIQUE:-0} ))
    if [[ "$needed" -le 0 ]]; then needed=1; fi
    printf "${ORANGE}${BOLD}pending${R}  ${GRAY}need %d more distinct frontier size(s)${R}" \
      "$needed"
  fi
}

# ── Live single-line display ──────────────────────────────────────────────────
print_live() {
  read_stats
  local bar
  bar=$(make_bar "$PCT")
  local elapsed=$(( $(date +%s) - START_TS ))
  local mm=$(( elapsed / 60 )) ss=$(( elapsed % 60 ))

  print_header
  printf "  ${GRAY}Started${R}  ${CYAN}$START_TIME${R}   ${GRAY}Elapsed${R}  ${CYAN}$(printf '%02d:%02d' $mm $ss)${R}   ${GRAY}Updated${R}  ${CYAN}$(date +%H:%M:%S)${R}\n\n"

  # Progress box
  hline
  rowf "N = ${YELLOW}${BOLD}${LAST_N:-?}${R} / 808"
  rowf "Lines = ${GREEN}${BOLD}${LINES:-?}${R}   Frontier = ${PURPLE}${BOLD}$(comma ${FRONTIER:-0})${R}"
  empty
  rowf "%b  ${YELLOW}${BOLD}%3d%%${R}" "$bar" "$PCT"
  bline

  printf "\n"
  # System + regression box
  hline
  rowf "RSS Memory   ${GREEN}${BOLD}$(human_kb ${RSS_KB:-0})${R}    Naive/task  ${YELLOW}${BOLD}${PER_TASK} KB${R}"
  rowf "Hardware     ${GRAY}${VCPUS} vCPUs   ${RAM_GB} GB RAM${R}"
  mline
  rowf "Regression   $(regression_status_line)"
  bline

  printf "\n${GRAY}  Waiting for N=808...   Press Ctrl-C to stop.${R}\n"
  printf   "${DIM}  Tip: running to N=851 (or until frontier ≥ 8,192 tasks) gives a${R}\n"
  printf   "${DIM}  more accurate asymptotic per-task measurement for the safety table.${R}\n"
}

# ── Full safety table ─────────────────────────────────────────────────────────
print_table() {
  read_stats
  local elapsed=$(( $(date +%s) - START_TS ))
  local mm=$(( elapsed / 60 )) ss=$(( elapsed % 60 ))

  print_header
  printf "  ${GRAY}Started${R}  ${CYAN}$START_TIME${R}   ${GRAY}Elapsed${R}  ${CYAN}$(printf '%02d:%02d' $mm $ss)${R}   ${GRAY}Updated${R}  ${CYAN}$(date +%H:%M:%S)${R}\n\n"

  # Progress box
  local bar
  bar=$(make_bar "$PCT")
  hline
  rowf "N = ${YELLOW}${BOLD}${LAST_N:-?}${R} / 808   Lines = ${GREEN}${BOLD}${LINES:-?}${R}   Frontier = ${PURPLE}${BOLD}$(comma ${FRONTIER:-0})${R}"
  empty
  rowf "%b  ${YELLOW}${BOLD}%3d%%${R}" "$bar" "$PCT"
  bline

  printf "\n"
  # System + regression box
  hline
  rowf "RSS Memory   ${GREEN}${BOLD}$(human_kb ${RSS_KB:-0})${R}    Naive/task  ${YELLOW}${BOLD}${PER_TASK} KB${R}"
  rowf "Hardware     ${GRAY}${VCPUS} vCPUs   ${RAM_GB} GB RAM${R}"
  mline
  rowf "Regression   $(regression_status_line)"
  bline

  printf "\n"

  # ── Safety table ──────────────────────────────────────────────────────────
  # Decide which per-task value to base predictions on, and build the header.
  local table_rate table_note
  if [[ "$ASYM_READY" -eq 1 ]]; then
    table_rate="$ASYM_BETA_KB"
    table_note="asymptotic β = ${ASYM_BETA_KB} KB/task"
  else
    table_rate="$PER_TASK"
    table_note="naive ${PER_TASK} KB/task  (regression pending)"
  fi

  hline
  rowf "${BOLD}MULTIPLIER SAFETY TABLE${R}  ${GRAY}(%s)${R}" "$table_note"
  empty
  rowf "${GRAY}%-10s  %-12s   %-12s   %-22s${R}" \
       "Multiplier" "Tasks" "Est. RAM" "Status"
  rowf "${DIM}%-10s  %-12s   %-12s   %-22s${R}" \
       "──────────" "──────────" "──────────" "────────────────────"

  for m in 131072 32768 16384 8192 2048 512 128 16 4; do
    local tasks=$(( VCPUS * m ))
    local est
    est=$(predict_rss_gb "$tasks")
    local icon col hint
    icon=$(safety_icon "$est" "$RAM_GB")
    col=$(safety_color  "$est" "$RAM_GB")
    if (( $(echo "$est < $RAM_GB" | bc -l) )); then
      hint="(fits in ${RAM_GB} GB)"
    else
      hint="(needs ${est} GB)"
    fi
    rowf "${GRAY}%8sU  →${R}  ${CYAN}%10s${R}  ${GRAY}→${R}  ${col}%8s GB${R}   ${col}%s${R} ${DIM}%s${R}" \
      "$(comma $m)" "$(comma $tasks)" "$est" "$icon" "$hint"
  done

  empty

  # ── Regression-quality callout ────────────────────────────────────────────
  if [[ "$ASYM_READY" -eq 1 ]]; then
    mline
    # Show the two end-points of the fit to give the user an intuition
    local f_small=$(( VCPUS * 512   ))
    local f_large=$(( VCPUS * 8192  ))
    local rss_small rss_large naive_small naive_large
    rss_small=$(predict_rss_gb "$f_small")
    rss_large=$(predict_rss_gb "$f_large")
    naive_small=$(echo "scale=2; $f_small * $PER_TASK / 1048576" | bc)
    naive_large=$(echo "scale=2; $f_large * $PER_TASK / 1048576" | bc)
    rowf "${GRAY}Model   per_task(F) = a + b/F  →  %d points, a≈%s KB/task, b≈%s MB${R}" \
      "$ASYM_UNIQUE" "$ASYM_BETA_KB" "$ASYM_ALPHA_MB"
    rowf "${GRAY}Predict RSS = a·F + b  (hyperbolic decay; accurate at all frontier sizes)${R}"
    rowf "${GRAY}512U  (%-7s tasks):  naive %5s GB  →  predicted %5s GB${R}" \
      "$(comma $f_small)" "$naive_small" "$rss_small"
    rowf "${GRAY}8192U (%-7s tasks):  naive %5s GB  →  predicted %5s GB${R}" \
      "$(comma $f_large)" "$naive_large" "$rss_large"
    rowf "${ORANGE}NOTE  Naive estimates are inflated (include shrinking per-task overhead).${R}"
    rowf "${GREEN}      Predicted estimates use asymptotic slope a (monotone, stable).${R}"
  else
    mline
    rowf "${ORANGE}NOTE  Regression not yet available — estimates use naive RSS/frontier.${R}"
    rowf "${ORANGE}      They are too conservative. Wait for a second distinct frontier${R}"
    rowf "${ORANGE}      size to appear in the log. The next increase occurs at N=851${R}"
    rowf "${ORANGE}      (cost 121, multiplier 512U → frontier = vcpus × 512 tasks).${R}"
  fi

  bline
  printf "\n  ${GRAY}Next refresh in ${CYAN}5s${GRAY}.   Press Ctrl-C to stop.${R}\n"
}

# =============================================================================
#  MAIN
# =============================================================================
cp "$SOURCE_PATH" ~/primecover.cpp
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions \
  -o primecover primecover.cpp || { echo "Build failed."; exit 1; }

(
  ./primecover > primecover.log 2>&1 &
  SOLVER_PID=$!
  START_TS=$(date +%s)
  START_TIME=$(date +%H:%M:%S)

  # Clean up sample file and restore terminal on exit
  trap 'rm -f "$SAMPLE_FILE"; kill $SOLVER_PID 2>/dev/null; tput cnorm; printf "\n\033[0m"; exit' INT

  tput civis   # hide cursor during live mode

  MODE="LIVE"

  while true; do
    if ! kill -0 "$SOLVER_PID" 2>/dev/null; then
      printf "\n${GREEN}${BOLD}  Solver finished.${R}\n\n"
      tput cnorm
      break
    fi

    read_stats   # populate globals + append sample + run regression

    if [[ "$MODE" = "LIVE" ]]; then
      if [[ -n "$LAST_N" && "$LAST_N" -ge 808 ]]; then
        MODE="TABLE"
        # Fall through immediately to print first table
      else
        print_live
        sleep 2
        continue
      fi
    fi

    # TABLE mode
    print_table
    sleep 5
  done

  tput cnorm
  rm -f "$SAMPLE_FILE"
  kill "$SOLVER_PID" 2>/dev/null
  wait "$SOLVER_PID" 2>/dev/null
)