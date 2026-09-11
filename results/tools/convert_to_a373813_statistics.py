# =============================================================================
# convert_to_a373813_statistics.py
# =============================================================================
# SPDX-License-Identifier: CC0-1.0
# This file is dedicated to the public domain. See ../LICENSE.
#
# Converts a373813_all_lines.txt into the clean a373813_statistics.txt format.
# (converts primecover2048_line_coordinates.cpp output to default output)
#
# a373813_all_lines.txt is a detailed log where each entry describes one prime
# number, containing a header line, a Stats block, and a Cover block. This
# script strips away the headers and cover sets, and extracts only the Stats
# content — one line per prime — producing a clean, readable statistics file.
#
# WR was ran on primecover2048_line_coordinates.cpp and converted to default output.
#
# Example output:
#   N=1 prime=2 lines=1 time=0.000396s mode=R active=0 ...
#   N=2 prime=3 lines=1 time=0.000224s mode=R active=0 ...
#   N=3 prime=5 lines=2 time=0.000177s mode=R active=0 ...
#
# -----------------------------------------------------------------------------
# USAGE (single command — update the path to match your setup):
#
#   cp /mnt/c/Math/prime-line-cover/results/tools/convert_to_a373813_statistics.py . && python3 convert_to_a373813_statistics.py /mnt/c/Math/prime-line-cover/results/a373813_all_lines.txt > /mnt/c/Math/prime-line-cover/results/a373813_statistics.txt
# =============================================================================

import re
import sys

def convert(input_path):
    with open(input_path, "r", encoding="utf-8") as f:
        content = f.read()

    # Split into individual prime blocks on the ── separator
    blocks = re.split(r'─{2,}', content)

    for block in blocks:
        # Extract everything between "Stats:" and "Cover:"
        stats_match = re.search(r'Stats:\s*(.*?)\s*Cover:', block, re.DOTALL)
        if stats_match:
            # Collapse any internal whitespace/newlines to a single space
            stats_line = re.sub(r'\s+', ' ', stats_match.group(1)).strip()
            print(stats_line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 convert_to_a373813_statistics.py a373813_all_lines.txt", file=sys.stderr)
        sys.exit(1)
    convert(sys.argv[1])