# =============================================================================
# convert_to_b373813.py
# =============================================================================
# SPDX-License-Identifier: CC0-1.0
# This file is dedicated to the public domain. See ../LICENSE.
#
# Converts a373813_all_lines.txt into the compact b373813.txt format.
#
# a373813_all_lines.txt is a detailed log where each entry describes one prime
# number: how many lines (intervals) were needed to cover all integers up to
# that prime. Each entry contains a Stats block with fields like N= and lines=.
#
# WR was ran on primecover2048_line_coordinates.cpp and converted to the b file.
#
# This script reads through that file, extracts the index (N) and line count
# (lines) from every Stats block, and outputs one row per prime:
#
#   <N> <lines>
#
# Example output:
#   1 1
#   2 1
#   3 2
#
# This produces the clean two-column sequence used for OEIS entry A373813.
#
# -----------------------------------------------------------------------------
# USAGE (single command — update the path to match your setup):
#
#   cp /mnt/c/Math/prime-line-cover/results/tools/convert_to_b373813.py . && python3 convert_to_b373813.py /mnt/c/Math/prime-line-cover/results/a373813_all_lines.txt > /mnt/c/Math/prime-line-cover/results/b373813.txt
# =============================================================================

import re
import sys

def convert(input_path):
    with open(input_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Split into individual prime blocks on the ── separator
    blocks = re.split(r'─{2,}', content)

    for block in blocks:
        n_match = re.search(r'\bN=(\d+)\b', block)
        lines_match = re.search(r'\blines=(\d+)\b', block)
        if n_match and lines_match:
            print(f"{n_match.group(1)} {lines_match.group(1)}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 convert_to_b373813.py a373813_all_lines.txt", file=sys.stderr)
        sys.exit(1)
    convert(sys.argv[1])