#!/usr/bin/env python3

import argparse
import sys
from collections import Counter

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Create read-level QC plots from SeqKit fx2tab output."
    )

    parser.add_argument(
        "--figure",
        required=True,
        help="Output PNG figure.",
    )

    parser.add_argument(
        "--summary",
        required=True,
        help="Output TSV summary.",
    )

    parser.add_argument(
        "--genome-size",
        type=int,
        required=True,
        help="Estimated haploid genome size in base pairs.",
    )

    return parser.parse_args()


args = parse_arguments()

length_counts = Counter()
quality_counts = Counter()
gc_counts = Counter()

number_of_reads = 0
total_bases = 0
quality_sum = 0.0
gc_sum = 0.0


for line in sys.stdin:
    line = line.strip()

    if not line or line.startswith("#"):
        continue

    columns = line.split("\t")

    if len(columns) < 4:
        continue

    read_length = int(columns[1])
    gc_percent = float(columns[2])
    average_quality = float(columns[3])

    # Exact length counts let us calculate an exact read N50.
    length_counts[read_length] += 1

    # Quality and GC are grouped into bins for compact plotting.
    quality_bin = int(average_quality)
    gc_bin = int(gc_percent)

    quality_counts[quality_bin] += 1
    gc_counts[gc_bin] += 1

    number_of_reads += 1
    total_bases += read_length
    quality_sum += average_quality
    gc_sum += gc_percent


if number_of_reads == 0:
    raise ValueError("No valid read records were received.")


mean_length = total_bases / number_of_reads
mean_read_quality = quality_sum / number_of_reads
mean_read_gc = gc_sum / number_of_reads
estimated_coverage = total_bases / args.genome_size


# Calculate read N50 from the exact length counts.
half_of_bases = total_bases / 2
cumulative_bases = 0
read_n50 = 0

for length in sorted(length_counts, reverse=True):
    cumulative_bases += length * length_counts[length]

    if cumulative_bases >= half_of_bases:
        read_n50 = length
        break


# Prepare length distribution and cumulative retained yield.
lengths = np.array(sorted(length_counts))
counts = np.array([length_counts[length] for length in lengths])
bases_at_length = lengths * counts

cumulative_bases_from_longest = np.cumsum(bases_at_length[::-1])[::-1]
retained_yield_percent = cumulative_bases_from_longest / total_bases * 100


quality_bins = np.array(sorted(quality_counts))
quality_values = np.array(
    [quality_counts[value] for value in quality_bins]
)

gc_bins = np.array(sorted(gc_counts))
gc_values = np.array(
    [gc_counts[value] for value in gc_bins]
)


# Write the compact summary table.
with open(args.summary, "w", encoding="utf-8") as output:
    output.write("metric\tvalue\n")
    output.write(f"number_of_reads\t{number_of_reads}\n")
    output.write(f"total_bases\t{total_bases}\n")
    output.write(f"mean_read_length_bp\t{mean_length:.2f}\n")
    output.write(f"read_n50_bp\t{read_n50}\n")
    output.write(f"mean_read_quality\t{mean_read_quality:.2f}\n")
    output.write(f"mean_read_gc_percent\t{mean_read_gc:.2f}\n")
    output.write(f"estimated_coverage_x\t{estimated_coverage:.2f}\n")


# Create a four-panel figure.
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

axes[0, 0].bar(
    lengths / 1000,
    counts,
    width=0.2,
    color="#4472C4",
)
axes[0, 0].axvline(
    read_n50 / 1000,
    color="red",
    linestyle="--",
    label=f"N50 = {read_n50 / 1000:.1f} kb",
)
axes[0, 0].set_title("Read-length distribution")
axes[0, 0].set_xlabel("Read length (kb)")
axes[0, 0].set_ylabel("Number of reads")
axes[0, 0].legend()

axes[0, 1].bar(
    quality_bins,
    quality_values,
    width=0.9,
    color="#70AD47",
)
axes[0, 1].set_title("Per-read average quality")
axes[0, 1].set_xlabel("Mean read quality (Phred)")
axes[0, 1].set_ylabel("Number of reads")

axes[1, 0].plot(
    lengths / 1000,
    retained_yield_percent,
    color="#ED7D31",
    linewidth=2,
)
axes[1, 0].axhline(
    60 / estimated_coverage * 100,
    color="red",
    linestyle="--",
    label="Yield required for approximately 60×",
)
axes[1, 0].set_title("Yield retained after a minimum-length cutoff")
axes[1, 0].set_xlabel("Minimum read length (kb)")
axes[1, 0].set_ylabel("Total bases retained (%)")
axes[1, 0].set_ylim(0, 105)
axes[1, 0].legend()

axes[1, 1].bar(
    gc_bins,
    gc_values,
    width=0.9,
    color="#A5A5A5",
)
axes[1, 1].set_title("Per-read GC distribution")
axes[1, 1].set_xlabel("GC (%)")
axes[1, 1].set_ylabel("Number of reads")

fig.suptitle(
    "PacBio HiFi read quality control\n"
    f"{number_of_reads:,} reads | "
    f"{total_bases / 1e9:.1f} Gb | "
    f"{estimated_coverage:.1f}× estimated coverage",
    fontsize=15,
)

fig.tight_layout()
fig.savefig(args.figure, dpi=200, bbox_inches="tight")
