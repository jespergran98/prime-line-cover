# =============================================================================
# convert_to_a373813_statistics.py
# =============================================================================
# Converts A373813_ALL_LINES.txt into the clean A373813_STATISTICS.txt format.
# (converts primecover1024_line_coordinates.cpp output to default output)
#
# A373813_ALL_LINES.txt is a detailed log where each entry describes one prime
# number, containing a header line, a Stats block, and a Cover block. This
# script strips away the headers and cover sets, and extracts only the Stats
# content — one line per prime — producing a clean, readable statistics file.
#
# Example output:
#   N=1 prime=2 lines=1 time=0.000172s mode=R active=0 ...
#   N=2 prime=3 lines=1 time=0.000092s mode=R active=0 ...
#   N=3 prime=5 lines=2 time=0.000068s mode=R active=0 ...
#
# -----------------------------------------------------------------------------
# USAGE (single command — update the path to match your setup):
#
#   cp /mnt/c/Math/prime-line-cover/results/convert_to_a373813_statistics.py . && python3 convert_to_a373813_statistics.py /mnt/c/Math/prime-line-cover/results/A373813_ALL_LINES.txt > /mnt/c/Math/prime-line-cover/results/A373813_STATISTICS.txt
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
        print("Usage: python3 convert_to_a373813_statistics.py A373813_ALL_LINES.txt", file=sys.stderr)
        sys.exit(1)
    convert(sys.argv[1])