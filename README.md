# Prime Line Cover – Exact Minimum Line Cover Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Live Demo](https://img.shields.io/badge/Live%20Demo-View-blue)](https://prime-line-cover.vercel.app)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20096556.svg)](https://doi.org/10.5281/zenodo.20096556)

**[Quick Start](#quick-start) · [File Descriptions](#file-descriptions) · [Interactive Demo](#interactive-demo) · [Paper](#paper) · [Results](#results) · [Google Cloud](#google-cloud) · [Citation](#citation)**

---

## Problem

Imagine you write down the first few prime numbers:

| Index | Prime |
|-------|-------|
| 1     | 2     |
| 2     | 3     |
| 3     | 5     |
| 4     | 7     |
| 5     | 11    |
| 6     | 13    |
| …     | …     |

Now plot these as points on a piece of graph paper:  
point 1 is at (1, 2), point 2 is at (2, 3), point 3 at (3, 5), point 4 at (4, 7), and so on.

**The goal:** Draw as few straight lines as possible so that every point lies on at least one line.

- A straight line can pass through **one point, two points, or many points** – whatever you need.
- You may use lines in any direction – vertical, horizontal, or slanted.
- The only restriction is that each line must be straight (no curves).
- Lines may overlap, and multiple lines can cross the same prime point.

**Example:** For the first 5 points (N = 5), you can cover all of them with just **2 lines**:

- One line passes through (1, 2) and (5, 11).
- Another line passes through (2, 3), (3, 5), and (4, 7) – three points on one straight line.

That's the best possible – you cannot cover all five points with only 1 line. So the answer for N = 5 is **2**. This is the last time two lines ever suffice: from N = 6 onward, no cover with fewer than three lines exists.

## Introduction to the Problem

If you're new to the problem, this Numberphile video by Brady Haran — featuring Neil Sloane, founder of the OEIS — is a great place to start:

<a href="https://www.youtube.com/watch?v=VFoIPlUalRY">
  <img src="https://img.youtube.com/vi/VFoIPlUalRY/maxresdefault.jpg" alt="Awkward Primes – Numberphile" width="400">
</a>

In the video, Brady Haran and Neil Sloane (founder of the OEIS) call the primes that force a new line **awkward primes** ([OEIS A393445](https://oeis.org/A393445)) — with a particularly stubborn one earning the nickname the **"party-pooper prime"**.

---

## Performance

Despite the enormous number of possible line combinations (which grows exponentially with N), this solver computes the optimal cover for the first 800 prime points in **less than 60 seconds** on a 10-year-old laptop.

On a Google Cloud `c4d-highcpu-8` instance (8 vCPUs, 15 GB RAM), it reaches the boundary of the old record at N = 861 in **22 minutes and 47 seconds** – a 744-fold speedup, obliterating that prior certified record, which required **282 hours, 26 minutes, and 31.5 seconds** using a general-purpose mixed-integer programming (MIP) solver.

For larger N the cost concentrates in a handful of instances: N = 1811 alone took **240 hours and 10 minutes** on the same hardware, 52% of the whole run – but the incremental sweep is so efficient that most N are solved in milliseconds. The current certified world record stands at **N = 1814**, where f(1814) = 229 and p_1814 = 15551, also computed on that machine. That sweep ran **461 hours, 35 minutes, and 22 seconds** – 19.2 days – and stopped only because the $300 free-trial credit ran out; the first 1024 indices, the entire previous record, account for just 47.5 of those hours.

A concrete illustration of that efficiency: the longest known plateau, f = 145 holding for 237 consecutive indices (N = 1079–1315), contains 233 primes that silently fall onto an existing optimal line — all 233 certified in a combined **2.74 seconds** — and its closer is p_1316 = 10831. The earlier plateau at f = 69 (N = 464–575) is now only the second longest; its closer p_576 = 4211 keeps the nickname the "party-pooper prime". The full breakdown across N = 1–1814: **1,141 witness hits** (constant-time certification), **362 root-only closures** (no branching), **311 full branch-and-bound searches** — 62.9%, 20.0% and 17.1% of the indices, with the branch-and-bound searches accounting for 99.87% of the wall-clock. Mode assignment is build-dependent and is not part of the optimality certificate: over N = 1–1024 this run recorded 613 / 235 / 176, differing from the published 1024 run at seven indices, while f(N) agrees at all 1024 indices. On all 262 non-witness steps in the new N = 1025–1814 block, the root gap ub − lb was exactly 1 before the first branch.

**What this repository does:** It finds the **exact minimum number of lines** needed for any N up to 1814. The full certified results for N = 1 through 1814 are included. The sequence of these minimum numbers is called [A373813](https://oeis.org/A373813) in the On-Line Encyclopedia of Integer Sequences (OEIS); the indices at which f(N) increases — the awkward primes — are catalogued as [A393445](https://oeis.org/A393445).

## Solver Features

Three versions are provided:

- **`primecover1024.cpp`** – core solver, outputs compact statistics for each N (time, search nodes, bounds, etc.).
- **`primecover1024_line_coordinates.cpp`** – identical solver that additionally writes out the exact coordinates of every line in the optimal cover (for full reproducibility and visualisation).
- **`primecover2048_line_coordinates.cpp`** – the current record solver: the line-coordinates solver with the coverage bitmasks widened from 1024 to 2048 bits (114 of 3,432 lines touched). The algorithm is byte-for-byte unchanged; there is no stats-only 2048 variant.

All three implement the exact algorithm described in the accompanying paper [`pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf`](https://prime-line-cover.vercel.app/pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf); the mask widening and the N = 1814 run are documented in the supplement [`pdf_extending_the_minimum_line_cover_to_n_1814.pdf`](https://prime-line-cover.vercel.app/pdf_extending_the_minimum_line_cover_to_n_1814.pdf). Key features:

- **Integer arithmetic throughout** – collinearity is checked via one integer multiply-compare; coverage counting is popcount over word-wise AND. No floating point anywhere on the hot path.
- **Heavy-line enumeration** – all lines containing at least three points of the final horizon are pre-computed once.
- **Bitmask representation** – supports up to N = 2048 points using 32x64-bit words in `primecover2048_line_coordinates.cpp` (the `1024` solvers use 16 words).
- **Exclusive Dependency Rule** – unconditional forcing of certain heavy lines (proved optimal).
- **Lagrangian relaxation** – provides tight lower bounds for pruning, with projected subgradient ascent and coordinate-descent polish.
- **Frontier decomposition** – the search tree is split into many independent tasks and processed in parallel over all CPU cores.
- **Incremental warm-start** – carries a witness cover, a warm heavy-line set, and a warm dual seed from one N to the next, making the sweep extremely efficient.

Each solver is optimised for the N in its filename: **N <= 1024** for the `1024` pair, **N <= 2048** for `primecover2048_line_coordinates.cpp`. The certified sweep stops at N = 1814 because the Google Cloud free-trial credit ran out, not because of any limit in the source — the 2048-bit ceiling is still 234 indices away and reaching it needs only more compute.

---

## Interactive Demo

**Live at: [https://prime-line-cover.vercel.app](https://prime-line-cover.vercel.app)**

`index.html` in the root of this repository is a complete single-file interactive demo. It is deployed as a static site on Vercel and requires no server – just open it in any modern browser. The JavaScript version is able to reach the 300th prime in 20 seconds (increase speed with the slider), but struggles beyond N = 700 (compared to the C++ version getting to N = 800 in 35-70 seconds depending on hardware).

The demo runs a faithful **JavaScript port of the C++ solver** (using the same bitmask-based algorithm, heavy-line enumeration, and branch-and-bound logic) directly in a Web Worker so the UI never blocks. It then animates the result step by step, incrementing N one prime at a time.

The JavaScript port is **single-threaded** and therefore slower than the native C++ solver for large N. It is well-suited for exploring and visualising the problem up to moderate N in real time, but is not intended for record-breaking runs. To compute the full certified sequence, use the C++ solver compiled locally.

---

## Paper

**[`pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf`](pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf)** — 32 pages. Still accurate about the architecture; only its numbers are superseded by the N = 1814 run.

**[`pdf_extending_the_minimum_line_cover_to_n_1814.pdf`](pdf_extending_the_minimum_line_cover_to_n_1814.pdf)** — 12 pages. The supplement covering the extension of the record to N = 1814 (f(1814) = 229, p = 15551): the 2048-bit mask widening, the 461 h 35 min sweep, and the 790 newly certified terms.

The main paper describes every component of the solver from mathematical foundations to implementation details. Here is exactly what each section covers:

**§1 Introduction** — States the problem, summarises the prior certified record (N = 861, 282 hours via MIP), lists the paper's concrete outcomes (f(1024) = 143 certified, 163 new OEIS terms, 20 new awkward primes, optimality proven for every reported value), sketches the NP-completeness of the general problem, and gives a section-by-section roadmap.

**§2 Problem Definition** — Defines prime points, heavy lines (lines through ≥ 3 prime points), the minimum cover number f(N), and the gain reformulation (maximising gain is equivalent to minimising line count). Also defines the full task structure representing a search node: the 1024-bit uncovered-point bitmask, residual coverage counts per heavy line, availability flags, accumulated gain, and two auxiliary point-indexed arrays that drive the Exclusive Dependency Rule.

**§3 Algorithm Overview** — Describes the two-phase architecture: (1) one-time heavy-line enumeration (~50 ms, producing 12,162 lines each stored as a 1024-bit bitmask, plus incidence and activation lookup structures), and (2) the incremental sweep. Covers the exact search skeleton, how the parallel frontier is built and expanded, and how independent subproblems are dispatched to worker threads via a shared atomic task queue. Includes Algorithms 1–3 (EnumerateHeavyLines, ApplyExclusiveDependencyRule, BuildFrontier) and the algorithm flow diagram.

**§4 Incremental Warm-Start and Witness Propagation** — Details the five items of carried state (certified cost, active heavy-line set, warm primal solution, warm dual seed, and witness cover), the three execution modes (W: witness hit certifies f(N) = f(N−1) without search; R: root bounds close the gap without branching; D: full parallel branch-and-bound), witness construction, and warm seed write-back.

**§5 Lower Bound** — Derives all bounding machinery: the cover-inequality lower bound ⌈|S|/2⌉ (from a greedily constructed 2-cover-free set S), the monotonicity floor f(N−1), the Lagrangian gain bound (an upper bound on attainable gain, computed via projected subgradient descent on per-point dual variables y ∈ [0,1]^N with a coordinate-descent polish pass), and branch-specific child-cap adjustments that derive tighter bounds for each branch direction individually.

**§6 Correctness and Optimality** — Proves four results: the Monotonicity Lemma (f(N+1) ≥ f(N)), the Exclusive Dependency Rule (if a productive heavy line has ≥ 3 points covered by no other productive line, forcing it never worsens the optimum), the Frontier Dominance Lemma (discarding nodes dominated in both uncovered-point set and blocked productive set preserves an optimal solution), and the Exact Solve Theorem, with mode-by-mode propositions (W, R, D each certify the correct f(N)) and a final certification corollary (every printed answer is proven optimal).

**§7 OEIS Sequence and Data** — Records the 163 newly certified terms of A373813 (N = 862–1024), the 20 new awkward primes ([A393445](https://oeis.org/A393445)), a staircase plot of f(N) for N = 1–1024, a plateau summary table (cover size, N-range, D/R/W row counts, wall-clock time per plateau), and a per-awkward-prime DFS diagnostics table (wall-clock time, active heavy lines, node count, Lagrangian iterations, nodes pruned by Lagrangian bound, strong-branching invocations, and maximum search depth).

**§8 Code Availability** — Documents the repository contents, the sweep log format and the meaning of every logged field, the frontier multiplier calibration ladder (including how the frontier diagnosis tool fits a RAM model to select the safe maximum multiplier per machine), and the precise architectural meaning of "1024" in the filename — and why extending the record requires widening the coverage masks to 2048 bits.

---

## Results

The repository contains result files in the `results/` folder, recording the computed cover costs (all for N = 1 to 1814):

- **`b373813.txt`** – simple two-column format: `N` and the optimal number of lines (easy for plotting), 1814 rows. Part of [OEIS A373813](https://oeis.org/A373813).
- **`a373813_statistics.txt`** – detailed statistics for each N (time, search nodes, bounds, etc.), one line per N for N = 1..1814. Derived from `a373813_all_lines.txt` by `results/tools/convert_to_a373813_statistics.py`.
- **`a373813_all_lines.txt`** – the full list of line coordinates (in `(x, y)` format) for each optimal cover, N = 1..1814 (~24 MB). The authoritative record, generated by `primecover2048_line_coordinates.cpp`.
- **`b393445.txt`** – simple two-column format: sequential index and prime value, listing only the awkward primes where the line count increases: 229 entries, ending at index 229, p = 15511. Part of [OEIS A393445](https://oeis.org/A393445).

The superseded N = 1024 record is preserved byte-for-byte in `results/archive_1024/`, under the same four filenames, so the two records diff directly.

Three Python converter scripts are also included in `results/tools/` to derive the above files from the raw solver output:

- **`convert_to_b373813.py`** – reads `a373813_all_lines.txt` and extracts the `N` and `lines` fields from every Stats block, producing the compact two-column `b373813.txt`.
- **`convert_to_a373813_statistics.py`** – reads `a373813_all_lines.txt` and extracts the full Stats line for every prime, producing the clean one-line-per-prime `a373813_statistics.txt`.
- **`convert_to_b393445.py`** – reads `a373813_all_lines.txt`, tracks the running maximum line count, and outputs only the primes where the line count strictly increases, producing `b393445.txt`.

---

## Quick Start

The solver uses C++23 features that only GCC 14 and newer support. Choose your operating system below.

⚠️ The [live demo](https://prime-line-cover.vercel.app) also contains this guide with one-click copy buttons and editable file paths — it may be easier to follow than the steps here.

### Windows (via WSL)

WSL (Windows Subsystem for Linux) lets you run a full Linux environment directly inside Windows — no virtual machine needed. It is the recommended way to build and run the solver on Windows.

**Step 1 — Install WSL**

Open **PowerShell as Administrator** (right-click the Start menu → "Windows PowerShell (Admin)") and run:

```powershell
wsl --install
```

This installs WSL 2 and the Ubuntu distribution in one go. Restart your PC when prompted, then open the **Ubuntu** app from the Start menu and wait for it to finish setting up your user account.

> Already have WSL and Ubuntu installed? Skip straight to Step 2.

**Step 2 — Install GCC 14**

The solver uses C++23 features that only GCC 14 and newer support. Ubuntu doesn't ship it by default, so these three commands add the right package repository and install it. You can paste all three at once — they will run in sequence:

```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt update
sudo apt install g++-14 -y
```

**Step 3 — Copy, compile, and run**

Download `primecover2048_line_coordinates.cpp` and save it to a location you'll remember — your Downloads folder or Desktop works well. Then replace the path in the command below with where you saved it, and paste the whole thing into your Ubuntu terminal:

```bash
cp /mnt/c/Users/YourName/Downloads/primecover2048_line_coordinates.cpp ~/solver.cpp && \
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  solver.cpp -o solver && \
./solver
```

> **Tip:** Your Windows `C:\` drive is always accessible inside WSL at `/mnt/c/`. So `C:\Users\Alice\Downloads` becomes `/mnt/c/Users/Alice/Downloads`.

The solver starts printing results immediately — one line per N.

For the previous 1024-bit solvers, replace `primecover2048_line_coordinates.cpp` in the `cp` command above with `primecover1024.cpp` (statistics only) or `primecover1024_line_coordinates.cpp` — the rest of the command stays the same.

---

### Linux (Ubuntu / Debian · Arch / Manjaro · Fedora / RHEL)

**Step 1 — Open a terminal**

Press `Ctrl+Alt+T` or launch a terminal from your application menu. No extra layer needed — you're running Linux natively.

**Step 2 — Install GCC 14**

The solver uses C++23 features that only GCC 14 and newer support. Run the command for your distribution:

*Ubuntu / Debian* — adds the toolchain PPA and installs `g++-14` (paste all three at once):

```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt update
sudo apt install g++-14 -y
```

*Arch / Manjaro* — GCC on a current Arch install is already GCC 14:

```bash
sudo pacman -S gcc
```

*Fedora / RHEL* — GCC on Fedora 40+ is already GCC 14:

```bash
sudo dnf install gcc-c++ -y
```

**Step 3 — Copy, compile, and run**

Download `primecover2048_line_coordinates.cpp` and save it somewhere in your home directory — your Downloads folder works well. Then replace the path below with where you saved it and paste the whole command into your terminal:

*Ubuntu / Debian*

```bash
cp /home/username/Downloads/primecover2048_line_coordinates.cpp ~/solver.cpp && \
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  solver.cpp -o solver && \
./solver
```

*Arch / Manjaro* — the compiler binary is `g++`, not `g++-14`:

```bash
cp /home/username/Downloads/primecover2048_line_coordinates.cpp ~/solver.cpp && \
g++ -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  solver.cpp -o solver && \
./solver
```

*Fedora / RHEL* — likewise, use `g++`:

```bash
cp /home/username/Downloads/primecover2048_line_coordinates.cpp ~/solver.cpp && \
g++ -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  solver.cpp -o solver && \
./solver
```

> Replace `username` with your actual Linux username. Run `whoami` in your terminal if you're unsure what it is.

The solver starts printing results immediately — one line per N.

For the previous 1024-bit solvers, replace `primecover2048_line_coordinates.cpp` in the `cp` command above with `primecover1024.cpp` (statistics only) or `primecover1024_line_coordinates.cpp` — the rest of the command stays the same.

---

### macOS

**Step 1 — Install Homebrew**

If you don't have Homebrew yet, open **Terminal** (⌘ Space → "Terminal") and run:

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

> Already have Homebrew? Skip straight to Step 2.

**Step 2 — Install GCC 14**

Apple Clang does not yet support all C++23 features used by the solver — GCC 14 via Homebrew is required. Install it with:

```bash
brew install gcc@14
```

Homebrew installs `g++-14` into its bin directory and adds it to your PATH automatically.

**Step 3 — Copy, compile, and run**

Download `primecover2048_line_coordinates.cpp` and save it somewhere in your home folder — your Downloads folder works well. Then replace the path below with where you saved it and paste the whole command into Terminal:

```bash
cp /Users/username/Downloads/primecover2048_line_coordinates.cpp ~/solver.cpp && \
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  solver.cpp -o solver && \
./solver
```

> Replace `username` with your macOS username. Run `whoami` in Terminal if you're unsure what it is.

The solver starts printing results immediately — one line per N.

For the previous 1024-bit solvers, replace `primecover2048_line_coordinates.cpp` in the `cp` command above with `primecover1024.cpp` (statistics only) or `primecover1024_line_coordinates.cpp` — the rest of the command stays the same.

---

### Google Cloud <a name="google-cloud"></a>

Google Cloud offers a **$300 free trial (90 days)** for new accounts — the credit that funded the full N = 1–1814 record sweep, 461 hours and 35 minutes of continuous solving. The sweep stopped when the credit ran out, not because the solver reached a limit. This guide uses `primecover2048_line_coordinates.cpp` — the solver that set the N = 1814 record, and the only 2048-bit build — which outputs the exact line coordinates of each optimal cover alongside the standard statistics.

> **Recommended instance:** `c4d-highcpu-8` — 8 vCPUs, 15 GB RAM, AMD Turin — the instance used to set the N = 1814 world record, on the free trial.
When creating the VM, set *Series* to **C4D**, *Machine type* to **c4d-highcpu-8**, and *Minimum CPU platform* to **AMD Turin**.

**Step 1 — Install GCC 14**

Connect via SSH (Compute Engine → VM Instances → SSH). GCC 14 is not in the default Debian repos — pull it from trixie:

```bash
echo "deb https://deb.debian.org/debian trixie main" | sudo tee /etc/apt/sources.list.d/trixie.list
printf "Package: *\nPin: release n=trixie\nPin-Priority: 100\n" | sudo tee /etc/apt/preferences.d/trixie
sudo apt update
sudo apt install -t trixie libc-bin -y
sudo dpkg --configure -a
sudo apt install -t trixie g++-14 -y
```

**Step 2 — Upload and compile**

Use the cloud icon (top-right of the SSH window) → **Upload File** → select `primecover2048_line_coordinates.cpp`. Then compile:

```bash
g++-14 -std=c++23 -O3 -march=znver5 -pthread -fno-exceptions -fno-rtti \
  primecover2048_line_coordinates.cpp -o solver
```

**Step 3 — Launch**

```bash
nohup stdbuf -oL ./solver > output.log 2>/dev/null &
```

The solver runs fully detached. You can close the browser tab or shut down your PC at any time — it keeps running on Google's infrastructure.

**Monitoring**

```bash
pgrep -a solver                                               # confirm it is running
tail -n 100 output.log                                        # check recent output
tail -f output.log                                            # watch live (Ctrl+C to stop watching)
grep "Stats:" output.log | awk '{print $2, $4}'              # N and line count over time
```

To download results: cloud icon → **Download File** → `output.log`.

---

## Understanding the Output

Each line of output corresponds to one value of N and looks like this:

```
N=1015 prime=8081 lines=143 time=377.238s mode=D active=12012 nodes=14300557 frontier=16384 ub0=143 lb_floor=142 lb_cov=91 lb=142 gap=1 ub_raw=148 lb_raw=130 gap_raw=18 root_ub_frac=0.331 warm_size=130 greedy_heavy=125 prod=12012 pinc=40315 forced_root=0 forced=18750 lag_iters=363411620 polish_sweeps=2889110 lag_prune=6636420 strong_branch=14 depth=96
```

### Primary fields

| Field | Description |
|-------|-------------|
| `N` | The index of the current problem instance: the solver is covering the first N prime points, where the k-th point is at coordinates (k, p_k). |
| `prime` | The numerical value of the N-th prime p_N — the y-coordinate of the most recently added point. |
| `lines` | The certified exact optimal answer f(N): the minimum number of straight lines required to cover all N prime points. |
| `time` | The wall-clock time elapsed to fully solve this instance. |
| `mode` | How optimality was certified — see the **Modes** section below. |
| `active` | The number of active heavy lines at step N, where a heavy line passes through at least 3 prime points and becomes active once its third-smallest point index has been reached. |

### Modes

| Mode | Name | Description |
|------|------|-------------|
| `W` | Witness hit | The new point (N, p_N) already lies on a line in the witness cover constructed at step N−1 (a cover deliberately built with a look‑ahead over future points). Hence f(N) = f(N−1) is certified in constant time with no further search. |
| `R` | Root closed | The lower and upper bounds converged at the root node, certifying optimality without any branch‑and‑bound search (though root‑level bound calculations still run). |
| `D` | Depth‑first search | A full branch‑and‑bound search was required to certify the answer. |

### Bound and gap fields

The solver works by maintaining a lower bound (best proven minimum) and an upper bound (best known cover) and closing the gap between them.

| Field | Description |
|-------|-------------|
| `ub0` | The initial incumbent (upper bound on f(N)) established before any search: the minimum of N, f(N−1)+1, the greedy heuristic cost, and f(N−1) when a witness hit is detected. |
| `ub_raw` | The raw upper bound from the greedy heuristic alone: the number of heavy lines chosen greedily plus ⌈residual/2⌉ for the remaining uncovered points. |
| `lb_floor` | The lower-bound floor inherited from the previous step — equal to f(N−1) when a warm-start solution is available, since no optimal cover can be cheaper than the prior optimum, and zero otherwise. |
| `lb_cov` | The cover‑inequality lower bound: the size of a **large (maximal)** subset S of the **N prime points** such that no heavy line covers more than 2 members of S. The greedy removal heuristic yields a valid lower bound of ⌈|S|/2⌉ (though S may not be the maximum‑cardinality set). |
| `lb_raw` | The raw Lagrangian lower bound at the root before any floor adjustment: the integer lower bound derived directly from the root Lagrangian gain upper bound via ⌈(N − ⌊gain_ub⌋) / 2⌉. |
| `lb` | The final certified lower bound on f(N): the maximum of lb_raw, lb_floor, and lb_cov. |
| `gap` | The optimality gap at certification: ub0 − lb. A gap of zero means optimality was certified **without branch‑and‑bound search** (i.e., at the root, via witness hit or bound convergence). |
| `gap_raw` | The raw gap ub_raw − lb_raw, reflecting the spread between the greedy upper bound and the unadjusted Lagrangian lower bound before any warm‑start floor is applied. |
| `root_ub_frac` | The fractional part of the root Lagrangian gain upper bound (gain_ub − ⌊gain_ub⌋). Values close to 1.0 indicate that lb_raw is close to rounding up by one, meaning a slightly tighter dual solution would improve the bound. |

### Search statistics

These fields are non‑zero only in mode D (full branch‑and‑bound).

| Field | Description |
|-------|-------------|
| `nodes` | Total DFS branch‑and‑bound nodes explored across all parallel worker threads. |
| `frontier` | Number of independent subtasks generated by root pre‑expansion and dispatched to the parallel thread pool; 1 means the single‑threaded fallback ran. |
| `lag_iters` | Total projected subgradient ascent iterations executed across every Lagrangian bound evaluation in the solve, summed over all threads and nodes. |
| `polish_sweeps` | Total coordinate‑descent polish sweeps applied to the Lagrangian dual solution; each sweep is a single pass tightening all dual variables, triggered only when the bound is within 0.5 of the pruning threshold. |
| `lag_prune` | Total DFS nodes pruned because the Lagrangian gain upper bound proved the node could not improve upon the incumbent, summed across all threads. **Erratum:** in the N = 1814 sweep this counter was still a 32-bit `int` and overflowed at the single hardest index, logging `lag_prune=-1423934502` at N = 1811; the true value there is 2,871,032,794 and the sweep total is 10,562,920,000 rather than the 6,267,952,704 the log sums to. The counter is write-only diagnostics — it never feeds a bound or a prune decision — so no certified f(N) is affected. Widened to `long long` in `primecover2048_line_coordinates.cpp`; the published logs are left exactly as the run produced them. |
| `strong_branch` | Total DFS nodes at which strong branching was invoked. This requires: at least 2 candidate lines in the shortlist, the node gap ≤ 1, at least 64 points remain uncovered, and the current depth ≤ 6. Up to 4 candidate lines are scored by simulating the include branch to select the variable that best prunes both child subtrees. |
| `depth` | Maximum DFS recursion depth reached across all parallel worker threads. |
| `forced` | Total heavy lines unconditionally forced into the solution by the Exclusive Dependency Rule across all nodes (root preprocessing pass + all DFS nodes) and all threads — a line is forced when it is the sole productive heavy line available to at least 3 of its active points. |
| `forced_root` | Heavy lines forced by the Exclusive Dependency Rule during the **initial root preprocessing pass only**, before the frontier is built or DFS begins. |

### Warm‑start and structure fields

These fields are valid for all modes (W, R, and D) and describe the problem structure and the knowledge carried in from the previous N.

| Field | Description |
|-------|-------------|
| `warm_size` | Number of heavy lines in the primal solution inherited from step N−1; zero means the solver started from a cold primal. |
| `greedy_heavy` | Number of heavy lines (each covering ≥ 3 points) selected by the greedy heuristic; ub_raw minus this gives the ⌈residual/2⌉ padding cost for points not reached by any chosen heavy line. |
| `prod` | Number of productive heavy lines at the root after forced‑line reduction: lines that are still available, cover at least 3 active points, and contribute a strictly positive gain. |
| `pinc` | Total productive incidence at the root: the sum of active‑point coverage counts over all productive lines, measuring how many (productive line, active point) pairs exist and thus how densely the productive lines overlap the active points. |

## Configuration Options

The solver's behaviour can be adjusted by modifying constants in the `config` namespace inside `primecover2048_line_coordinates.cpp`:

| Constant | Default | Description |
|----------|---------|-------------|
| `kBitCapacity`          | 2048 | Hard limit on N (bitset size). `kBitWords` is derived from it (2048 / 64 = 32). Do not increase beyond 2048 without widening the bitmask type. |
| `kStartN`               | 1    | Starting N for the sweep. Increase to resume mid-sequence. |
| `kExecutionLimit`       | 2048 | Maximum N to compute (capped by `kBitCapacity`). |
| `kPerNTimeLimitSeconds` | 0    | Per-instance time limit in seconds (0 = no limit). If exceeded, the solver prints `[stopped]` and exits. |

By default the solver runs from N = 1 up to N = 2048. To change that, edit `kStartN` and `kExecutionLimit` at the top of the source file. The published sweep stopped at N = 1814 because its compute budget ran out, not because of any limit in the source: the remaining 234 indices up to the 2048-bit ceiling need no source change, only more compute.

## Performance Tuning

The solver uses all available CPU cores (`std::thread::hardware_concurrency()`) automatically. The options below start to matter once the cost ladder reaches its larger rungs (N > ~850), and they matter most for an attempt past the current record at N = 1814.

### Frontier multiplier

The solver decomposes the search tree into independent tasks and distributes them across CPU cores. The number of tasks generated is:

```
frontier size  =  CPU threads  ×  frontier_multiplier
```

Larger multipliers increase parallelism but also memory usage. The multiplier is selected automatically via a cost ladder in `primecover2048_line_coordinates.cpp`:

```cpp
const unsigned frontier_multiplier =
//  current_best_cost() >= 143 ? 131072U  // rarely feasible
//: current_best_cost() >= 138 ?  32768U  // days–weeks
//: current_best_cost() >= 133 ?  16384U  // only if diagnostic confirms it is safe
//  current_best_cost() >= 128 ?   8192U  // needs ~30 GB with 256-byte masks; disabled for the 15 GB record run
    current_best_cost() >= 125 ?   2048U  // safe on c4d-highcpu-8 (15 GB) and on a personal i9 9900k
  : current_best_cost() >= 121 ?    512U  // safe on most machines
  : current_best_cost() >= 113 ?    128U  // safe
  : current_best_cost() >= 93  ?     16U  // sub-second
                               :      4U; // trivial
```

The defaults are conservative and safe. As shipped, the ladder tops out at `2048U` — 16,384 tasks on the 8-thread reference instance — and that is the ceiling the entire N = 1–1814 sweep ran under on 15 GB. **For an attempt past the current record at N = 1814 you should measure your machine's per-task RAM usage and unlock the highest multiplier that fits in your available memory.** Per-task cost varies significantly between machines — for the 1024-bit build, ~184 KB/task on the reference cloud instance and ~125 KB/task on some desktops — and the 2048-bit build roughly doubles that, since every heavy-line mask is 256 bytes instead of 128. Do not guess.

### Frontier diagnosis tool

`tools/primecover_frontier_diagnosis.sh` automates the measurement. It compiles and runs the solver, displays a live progress dashboard, and — once N reaches 808 — fits a hyperbolic decay model (`per_task(F) = a + b/F`) to the observed RSS samples and produces a safety table like this:

```
Multiplier   Tasks        Est. RAM     Status
──────────   ──────────   ──────────   ────────────────────
  131072U →  1,048,576  →  386.4 GB   ✗  UNSAFE
   32768U →    262,144  →   97.0 GB   ✗  UNSAFE
   16384U →    131,072  →   48.8 GB   ✗  UNSAFE
    8192U →     65,536  →   24.6 GB   ✗  UNSAFE
    2048U →     16,384  →    6.4 GB   ✓  SAFE
    ...
```

> The RAM estimates above are illustrative — your values will differ based on your hardware and thread count.

To use it:

1. Open `tools/primecover_frontier_diagnosis.sh` and set `SOURCE_PATH` on line 101 to the path of your `primecover2048_line_coordinates.cpp` file.
2. Run the script:
   ```bash
   bash tools/primecover_frontier_diagnosis.sh
   ```
3. When the safety table appears, find the highest multiplier marked **✓ SAFE**.
4. In `primecover2048_line_coordinates.cpp`, find the `frontier_multiplier` ladder and apply the edit shown below (example: unlocking `8192U`, the rung the record run had to leave disabled):

   ```cpp
   // Before:
   //  current_best_cost() >= 128 ?  8192U
       current_best_cost() >= 125 ?  2048U

   // After:
       current_best_cost() >= 128 ?  8192U
     : current_best_cost() >= 125 ?  2048U
   ```

   Remove the leading `//` from the line you are enabling (and its `:`, if it has one), then add a `:` before the line that was previously the first active condition. Recompile and run.

The diagnostic only needs to run once per machine. It is not needed for runs up to N ≈ 850 — the default settings are safe there — and as shipped the ladder tops out at `2048U`, the ceiling that carried the entire N = 1–1814 sweep on a 15 GB instance. Run it when you want to re-enable `8192U`, which in the 2048-bit build needs roughly 30 GB.

### CPU architecture flag

Compile with `-march=native` (already in the default command) to let GCC auto-detect your CPU. On **AMD Zen 5** machines you can be explicit for a small extra gain:

```bash
g++-14 -std=c++23 -O3 -march=znver5 -pthread -fno-exceptions -fno-rtti \
  solver.cpp -o solver
```

### WSL CPU priority (Windows only)

To give WSL higher CPU priority while the solver is running, open **PowerShell as Administrator** and execute:

```powershell
(Get-Process vmmem*).PriorityClass = 'High'
```

## File Descriptions

| File | Description |
|------|-------------|
| `primecover1024.cpp` | Solver source code for N <= 1024 (stats-only output). Used for the superseded N = 1024 record. |
| `primecover1024_line_coordinates.cpp` | Same solver, but outputs full line coordinates. Produced the superseded N = 1024 record. |
| `primecover2048_line_coordinates.cpp` | Current record solver: the same algorithm with 2048-bit masks (`kBitCapacity = 2048`, `kBitWords = 32`). Outputs full line coordinates; there is no stats-only variant. Produced the certified N = 1814 record. |
| `results/b373813.txt` | Two-column `N` and optimal lines (space-separated), N = 1..1814. Part of [OEIS A373813](https://oeis.org/A373813). |
| `results/a373813_statistics.txt` | Detailed statistics for each N (time, search nodes, bounds, etc.), N = 1..1814. Derived from `a373813_all_lines.txt` by `results/tools/convert_to_a373813_statistics.py`. |
| `results/a373813_all_lines.txt` | Full line-by-line coordinates of each optimal cover, N = 1..1814, from `primecover2048_line_coordinates.cpp` (~24 MB). The authoritative sweep log; the other three data files are derived from it. |
| `results/b393445.txt` | Two-column sequential index and prime value, listing only the awkward primes where the line count increases — 229 of them up to N = 1814. Part of [OEIS A393445](https://oeis.org/A393445). |
| `results/archive_1024/` | The superseded N = 1024 record, preserved byte-for-byte as published under the same four filenames, so the two records diff directly. |
| `results/LICENSE` | CC0 1.0 Universal dedication covering everything under `results/` — the data is public domain, unlike the MIT-licensed code. |
| `results/tools/convert_to_b373813.py` | Converter script: extracts `N` and `lines` from every Stats block in `a373813_all_lines.txt` to produce `b373813.txt`. |
| `results/tools/convert_to_a373813_statistics.py` | Converter script: extracts the full Stats line for every prime from `a373813_all_lines.txt` to produce `a373813_statistics.txt`. |
| `results/tools/convert_to_b393445.py` | Converter script: extracts only the primes where the line count strictly increases from `a373813_all_lines.txt` to produce `b393445.txt`. |
| `index.html` | Self-contained interactive demo (JavaScript port of the solver). Deploy as a static site or open locally — no build step required. See [Interactive Demo](#interactive-demo). |
| `pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf` | Mathematical paper describing the problem and the algorithm (32 pages). The architecture it documents is unchanged; its numbers describe the superseded N = 1024 record. |
| `pdf_extending_the_minimum_line_cover_to_n_1814.pdf` | Supplement paper (12 pages) documenting the extension of the certified record to N = 1814. |
| `tools/primecover_frontier_diagnosis.sh` | Frontier diagnosis script: compiles and runs the solver, displays a live dashboard, and produces a per-multiplier RAM safety table to help you unlock the highest safe frontier multiplier for your hardware. See [Performance Tuning](#performance-tuning). |
| `CITATION.cff` | Citation metadata in Citation File Format (CFF). Recognised automatically by GitHub, Zenodo, and Zotero. |
| `old_solvers/` | Earlier milestone variants of the solver, preserved for reproducibility. Provided as-is; not actively maintained. |

## Contributing

Bug reports and correctness challenges are welcome via [GitHub Issues](../../issues).

## About

Built by Jesper Gran Mikkelsen, an independent researcher in Norway with no prior background in operations research. The solver was developed over approximately 300 hours of AI-assisted iteration using Claude Code and ChatGPT (GPT-5.5 in Codex), prompted by watching the [Numberphile video on awkward primes](https://youtu.be/VFoIPlUalRY) in April 2026.

## Citation

If you use this code or the computed results in your work, please cite the accompanying papers — the original solver paper and the N = 1814 supplement that extends the record.

A [`CITATION.cff`](CITATION.cff) file is included in the repository — GitHub, Zenodo, and Zotero can read it automatically. For manual use, the BibTeX entry is:

```bibtex
@software{primecover2026,
  author    = {Jesper Gran Mikkelsen},
  title     = {prime-line-cover: Exact Minimum Line Cover Solver for Prime Points},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20096556},
  url       = {https://github.com/jespergran98/prime-line-cover}
}
```

## License

[MIT](LICENSE) – you are free to use, modify, and distribute the code, provided the original copyright and permission notice are included.