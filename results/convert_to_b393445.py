# =============================================================================
# convert_to_b393445.py
# =============================================================================
# SPDX-License-Identifier: CC0-1.0
# This file is dedicated to the public domain. See LICENSE in this directory.
#
# Converts A373813_ALL_LINES.txt into the compact B393445.txt format.
#
# A373813_ALL_LINES.txt is a detailed log where each entry describes one prime
# number: how many lines (intervals) were needed to cover all integers up to
# that prime. Each entry contains a Stats block with fields like prime= and lines=.
#
# An "awkward prime" is a prime where the line count increases compared to the
# previous prime — i.e., a new minimum number of lines is required for the
# first time at that prime.
#
# This script reads through that file, tracks the running maximum line count,
# and outputs one row for each prime where the line count is strictly greater
# than all previous primes:
#
#   <index> <prime>
#
# Example output:
#   1 2
#   2 5
#   3 13
#   4 29
#   5 59
#
# This produces the clean two-column sequence used for OEIS entry A393445.
#
# -----------------------------------------------------------------------------
# USAGE (single command — update the path to match your setup):
#
#   cp /mnt/c/Math/prime-line-cover/results/convert_to_b393445.py . && python3 convert_to_b393445.py /mnt/c/Math/prime-line-cover/results/A373813_ALL_LINES.txt > /mnt/c/Math/prime-line-cover/results/B393445.txt
# =============================================================================

import re
import sys

def convert(input_path):
    with open(input_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Split into individual prime blocks on the ── separator
    blocks = re.split(r'─{2,}', content)

    max_lines_seen = 0
    awkward_index = 0

    for block in blocks:
        prime_match = re.search(r'\bprime=(\d+)\b', block)
        lines_match = re.search(r'\blines=(\d+)\b', block)
        if prime_match and lines_match:
            prime = int(prime_match.group(1))
            lines = int(lines_match.group(1))
            if lines > max_lines_seen:
                max_lines_seen = lines
                awkward_index += 1
                print(f"{awkward_index} {prime}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 convert_to_b393445.py A373813_ALL_LINES.txt", file=sys.stderr)
        sys.exit(1)
    convert(sys.argv[1])