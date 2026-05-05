# Prime Line Cover – Exact Minimum Line Cover Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Live Demo](https://img.shields.io/badge/Live%20Demo-View-blue)](https://prime-line-cover.vercel.app)

**[Quick Start](#quick-start) · [File Descriptions](#file-descriptions) · [Interactive Demo](#interactive-demo) · [Results](#results) · [Citation](#citation)**

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

That's the best possible – you cannot cover all five points with only 1 line. So the answer for N = 5 is **2**.

## Introduction to the Problem

If you're new to the problem, this Numberphile video by Brady Haran — featuring Neil Sloane, founder of the OEIS — is a great place to start:

<a href="https://www.youtube.com/watch?v=VFoIPlUalRY">
  <img src="https://img.youtube.com/vi/VFoIPlUalRY/maxresdefault.jpg" alt="Awkward Primes – Numberphile" width="400">
</a>

---

## Performance

Despite the enormous number of possible line combinations (which grows exponentially with N), this solver computes the optimal cover for the first 800 prime points in **less than 60 seconds** on a 10-year-old laptop.

On a Google Cloud `c4d-highcpu-8` instance (8 vCPUs, 15 GB RAM), it reaches the previous world-record boundary at N = 861 in **just 22 minutes** – obliterating the prior certified record, which required **282 hours** using a general-purpose mixed-integer programming (MIP) solver.

For larger N up to 1024, the hardest instances take about an hour on the same hardware – but the incremental sweep is so efficient that most N are solved in milliseconds. The current certified world record stands at **N = 1024**, also computed on that machine.

**What this repository does:** It finds the **exact minimum number of lines** needed for any N up to 1024. The full certified results for N = 1 through 1024 are included. The sequence of these minimum numbers is called [A373813](https://oeis.org/A373813) in the On-Line Encyclopedia of Integer Sequences (OEIS).

## Solver Features

Two versions are provided:

- **`primecover1024.cpp`** – core solver, outputs compact statistics for each N (time, search nodes, bounds, etc.).
- **`primecover1024_line_coordinates.cpp`** – identical solver that additionally writes out the exact coordinates of every line in the optimal cover (for full reproducibility and visualisation).

Both implement the exact algorithm described in the accompanying paper [`pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf`](pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf). Key features:

- **Heavy-line enumeration** – all lines containing at least three points of the final horizon are pre-computed once.
- **Bitmask representation** – supports up to N = 1024 points using 16x64-bit words.
- **Exclusive Dependency Rule** – unconditional forcing of certain heavy lines (proved optimal).
- **Lagrangian relaxation** – provides tight upper bounds for pruning, with projected subgradient ascent and coordinate-descent polish.
- **Frontier decomposition** – the search tree is split into many independent tasks and processed in parallel over all CPU cores.
- **Incremental warm-start** – carries a witness cover, a warm heavy-line set, and a warm dual seed from one N to the next, making the sweep extremely efficient.

The solver is optimised for **N <= 1024** (hence the `1024` in the filename).

---

## Interactive Demo

**Live at: [https://prime-line-cover.vercel.app](https://prime-line-cover.vercel.app)**

`index.html` in the root of this repository is a complete single-file interactive demo. It is deployed as a static site on Vercel and requires no server – just open it in any modern browser. The JavaScript version is able to reach the 250th prime in 10 seconds (increase speed with the slider), but struggles beyond N = 280 (compared to the C++ version getting to N = 800 in 35-70 seconds depending on hardware).

The demo runs a faithful **JavaScript port of the C++ solver** (using the same bitmask-based algorithm, heavy-line enumeration, and branch-and-bound logic) directly in a Web Worker so the UI never blocks. It then animates the result step by step, incrementing N one prime at a time.

The JavaScript port is **single-threaded** and therefore slower than the native C++ solver for large N. It is well-suited for exploring and visualising the problem up to moderate N in real time, but is not intended for record-breaking runs. To compute the full certified sequence, use the C++ solver compiled locally.

---

## Results

The repository contains three files in the `/results` folder, recording the computed cover costs (all for N = 1 to 1024):

- **`B373813.txt`** – simple two-column format: `N` and the optimal number of lines (easy for plotting). Part of [OEIS A373813](https://oeis.org/A373813).
- **`A373813_STATISTICS.txt`** – detailed statistics for each N (time, search nodes, bounds, etc.). Generated by `primecover1024.cpp`.
- **`A373813_ALL_LINES.txt`** – the full list of line coordinates (in `(x, y)` format) for each optimal cover. Generated by `primecover1024_line_coordinates.cpp`.

The accompanying paper [`pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf`](pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf) describes the algorithm in full: the theoretical foundations of the Lagrangian relaxation bounds, the Exclusive Dependency Rule, and the incremental warm-start strategy. It includes per-instance statistics tables, a staircase plot of f(N), and a proof that every reported value is certified optimal. The role of **awkward primes** — indices at which the minimum line count increases — is discussed throughout.

---

## Quick Start

The solver uses C++23 features that only GCC 14 and newer support. Choose your operating system below.

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

Download `primecover1024.cpp` and save it to a location you'll remember — your Downloads folder or Desktop works well. Then replace the path in the command below with where you saved it, and paste the whole thing into your Ubuntu terminal:

```bash
cp /mnt/c/Users/YourName/Downloads/primecover1024.cpp ~/primecover.cpp && \
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  primecover.cpp -o primecover && \
./primecover
```

> **Tip:** Your Windows `C:\` drive is always accessible inside WSL at `/mnt/c/`. So `C:\Users\Alice\Downloads` becomes `/mnt/c/Users/Alice/Downloads`.

The solver starts printing results immediately — one line per N.

For the line-coordinates version, replace `primecover1024.cpp` with `primecover1024_line_coordinates.cpp` in the commands above.

---

### Linux (Ubuntu / Debian)

No extra setup needed — you're already running Linux natively. Open a terminal with `Ctrl+Alt+T` or from your application menu.

**Step 1 — Install GCC 14**

The solver uses C++23 features that only GCC 14 and newer support. Ubuntu and Debian don't ship it by default, so these three commands add the right package repository and install it. You can paste all three at once:

```bash
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt update
sudo apt install g++-14 -y
```

**Step 2 — Compile and run**

Download `primecover1024.cpp` and save it somewhere in your home directory — your Downloads folder works well. Then replace the path below with where you saved it and paste the whole command into your terminal:

```bash
cp /home/username/Downloads/primecover1024.cpp ~/primecover.cpp && \
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  primecover.cpp -o primecover && \
./primecover
```

> Replace `username` with your actual Linux username. Run `whoami` in your terminal if you're unsure what it is.

The solver starts printing results immediately — one line per N.

For the line-coordinates version, replace `primecover1024.cpp` with `primecover1024_line_coordinates.cpp` in the commands above.

---

### Linux (Arch / pacman)

Open a terminal and install GCC 14 with pacman:

**Step 1 — Install GCC 14**

```bash
sudo pacman -S gcc14
```

**Step 2 — Compile and run**

Download `primecover1024.cpp` and save it somewhere in your home directory. Then replace the path below with where you saved it and paste the whole command into your terminal:

```bash
cp /home/username/Downloads/primecover1024.cpp ~/primecover.cpp && \
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  primecover.cpp -o primecover && \
./primecover
```

> Replace `username` with your actual Linux username. Run `whoami` in your terminal if you're unsure what it is.

The solver starts printing results immediately — one line per N.

For the line-coordinates version, replace `primecover1024.cpp` with `primecover1024_line_coordinates.cpp` in the commands above.

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

Download `primecover1024.cpp` and save it somewhere in your home folder — your Downloads folder works well. Then replace the path below with where you saved it and paste the whole command into Terminal:

```bash
cp /Users/username/Downloads/primecover1024.cpp ~/primecover.cpp && \
g++-14 -std=c++23 -O3 -march=native -pthread -fno-exceptions -fno-rtti \
  primecover.cpp -o primecover && \
./primecover
```

> Replace `username` with your macOS username. Run `whoami` in Terminal if you're unsure what it is.

The solver starts printing results immediately — one line per N.

For the line-coordinates version, replace `primecover1024.cpp` with `primecover1024_line_coordinates.cpp` in the commands above.

> **Note:** For **AMD Zen 5** CPUs you can substitute `-march=native` with `-march=znver5` for a small extra performance gain.

---

## Understanding the Output

Each line of output looks like this:

```
N=862 prime=6689 lines=123 time=0.016016s mode=R active=9254 ...
```

- `lines` – the optimal cover cost f(N).
- `mode` – how the solution was obtained:
  - **W** = witness hit (inherited from previous prefix, constant-time)
  - **R** = root closed (bounds closed without DFS)
  - **D** = depth-first search (full branch-and-bound needed)
- `active` – number of active heavy lines at step N.
- Other fields (`ub0`, `lb_cov`, `lb`, `gap`, `nodes`, …) are detailed diagnostics – see the paper for their meaning.

For the line-coordinates version, each statistics line is followed by a `Cover:` block listing all lines in the optimal cover (each line as a list of `(x, y)` points).

## Configuration Options

The solver's behaviour can be adjusted by modifying constants in the `config` namespace inside `primecover1024.cpp`:

| Constant | Default | Description |
|----------|---------|-------------|
| `kBitCapacity`          | 1024 | Hard limit on N (bitset size). Do not increase beyond 1024 without changing the bitmask type. |
| `kStartN`               | 1    | Starting N for the sweep. Increase to resume mid-sequence. |
| `kExecutionLimit`       | 1024 | Maximum N to compute (capped by `kBitCapacity`). |
| `kPerNTimeLimitSeconds` | 0    | Per-instance time limit in seconds (0 = no limit). If exceeded, the solver prints `[stopped]` and exits. |

By default the solver runs from N = 1 up to N = 1024. To change that, edit `kStartN` and `kExecutionLimit` at the top of the source file.

## Performance Tuning

- The solver uses all CPU cores (`std::thread::hardware_concurrency()`).
- The **frontier multiplier** (see extensive comments in the source) controls how many parallel tasks are generated. This is tuned automatically based on the current best cost and your available RAM. For high-memory machines (>=30 GB), you can increase the top multiplier to `131072U` for better parallelism.
- To give WSL higher CPU priority on Windows, open **PowerShell as Administrator** while the solver is running and execute:
  ```powershell
  (Get-Process vmmem*).PriorityClass = 'High'
  ```

## File Descriptions

| File | Description |
|------|-------------|
| `primecover1024.cpp` | Main solver source code (stats-only output). |
| `primecover1024_line_coordinates.cpp` | Same solver, but outputs full line coordinates. |
| `results/B373813.txt` | Two-column `N` and optimal lines (space-separated), N = 1..1024. |
| `results/A373813_STATISTICS.txt` | Detailed statistics for each N (time, search nodes, bounds, etc.), from `primecover1024.cpp`. |
| `results/A373813_ALL_LINES.txt` | Full line-by-line coordinates of each optimal cover, from `primecover1024_line_coordinates.cpp`. |
| `index.html` | Self-contained interactive demo (JavaScript port of the solver). Deploy as a static site or open locally — no build step required. See [Interactive Demo](#interactive-demo). |
| `pdf_exact_solver_for_minimum_line_cover_of_prime_points.pdf` | Mathematical paper describing the problem and the algorithm. |
| `old_solvers/` | Earlier milestone variants of the solver, preserved for reproducibility. Provided as-is; not actively maintained. |

## Old Solvers

The `old_solvers/` folder contains earlier milestone variants of the solver used during development, preserved for reproducibility and to document the development history leading to the current version. They are provided as-is and are not actively maintained.

## Contributing

Bug reports and correctness challenges are welcome via [GitHub Issues](../../issues).

## Citation

If you use this code or the computed results in your work, please cite the accompanying paper:

```bibtex
@unpublished{primecover2026,
  author = {Jesper Gran Mikkelsen},
  title  = {An Incremental Exact Solver for the Minimum Line Cover of Prime Points},
  note   = {Unpublished manuscript},
  year   = {2026}
}
```

## License

[MIT](LICENSE) – you are free to use, modify, and distribute the code, provided the original copyright and permission notice are included.