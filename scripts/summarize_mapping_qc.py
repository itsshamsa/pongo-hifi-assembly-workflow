import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


input_file = "results/mapping/primary.coverage.tsv"
summary_file = "results/mapping/primary.mapping_summary.tsv"
flagged_file = "results/mapping/primary.flagged_contigs.tsv"
figure_file = "results/mapping/primary.mapping_qc.png"
parser = argparse.ArgumentParser(
    description="Summarize per-contig mapping coverage."
)

parser.add_argument("--coverage", required=True)
parser.add_argument("--summary", required=True)
parser.add_argument("--flagged", required=True)
parser.add_argument("--figure", required=True)

args = parser.parse_args()

input_file = args.coverage
summary_file = args.summary
flagged_file = args.flagged
figure_file = args.figure

# Read the tab-separated SAMtools coverage table
coverage = pd.read_csv(input_file, sep="\t")

# Rename the first column to make it easier to use
coverage = coverage.rename(columns={"#rname": "contig"})

# Calculate the length represented by each row
coverage["length_bp"] = (
    coverage["endpos"] - coverage["startpos"] + 1
)

total_bp = coverage["length_bp"].sum()
covered_bp = coverage["covbases"].sum()

genome_breadth = covered_bp / total_bp * 100

weighted_mean_depth = np.average(
    coverage["meandepth"],
    weights=coverage["length_bp"],
)

median_contig_depth = coverage["meandepth"].median()

low_depth_threshold = weighted_mean_depth * 0.5
high_depth_threshold = weighted_mean_depth * 2

flagged = coverage[
    (coverage["meandepth"] < low_depth_threshold)
    | (coverage["meandepth"] > high_depth_threshold)
    | (coverage["coverage"] < 99)
].copy()

flagged_bp = flagged["length_bp"].sum()
flagged_bp_percent = flagged_bp / total_bp * 100

low_depth_contigs = (
    coverage["meandepth"] < low_depth_threshold
).sum()

high_depth_contigs = (
    coverage["meandepth"] > high_depth_threshold
).sum()

depth_above_50_contigs = (
    coverage["meandepth"] > 50
).sum()

summary = pd.DataFrame(
    {
        "metric": [
            "contigs",
            "total_bp",
            "covered_bp",
            "genome_breadth_percent",
            "weighted_mean_depth",
            "median_contig_depth",
            "contigs_coverage_below_99_percent",
            "flagged_contigs",
        ],
        "value": [
            len(coverage),
            total_bp,
            covered_bp,
            genome_breadth,
            weighted_mean_depth,
            median_contig_depth,
            (coverage["coverage"] < 99).sum(),
            len(flagged),
        ],
    }
)

summary = pd.DataFrame(
    {
        "metric": [
            "contigs",
            "total_bp",
            "covered_bp",
            "genome_breadth_percent",
            "weighted_mean_depth",
            "median_contig_depth",
            "contigs_coverage_below_99_percent",
            "low_depth_contigs",
            "high_depth_contigs",
            "contigs_depth_above_50x",
            "flagged_contigs",
            "flagged_bp",
            "flagged_bp_percent",
        ],
        "value": [
            len(coverage),
            total_bp,
            covered_bp,
            genome_breadth,
            weighted_mean_depth,
            median_contig_depth,
            (coverage["coverage"] < 99).sum(),
            low_depth_contigs,
            high_depth_contigs,
            depth_above_50_contigs,
            len(flagged),
            flagged_bp,
            flagged_bp_percent,
        ],
    }
)

summary.to_csv(summary_file, sep="\t", index=False)

flagged.sort_values("meandepth").to_csv(
    flagged_file,
    sep="\t",
    index=False,
)


# Create the QC figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

scatter = axes[0].scatter(
    coverage["length_bp"] / 1_000_000,
    coverage["meandepth"],
    c=coverage["meanmapq"],
    cmap="viridis",
    alpha=0.7,
    s=20,
)

axes[0].axhline(
    weighted_mean_depth,
    color="red",
    linestyle="--",
    label=f"Weighted mean: {weighted_mean_depth:.2f}×",
)

axes[0].set_xscale("log")
axes[0].set_xlabel("Contig length (Mb, log scale)")
axes[0].set_ylabel("Mean depth")
axes[0].set_title("Coverage depth by contig")
axes[0].legend()

colorbar = fig.colorbar(scatter, ax=axes[0])
colorbar.set_label("Mean mapping quality")

axes[1].hist(
    coverage["meandepth"],
    bins=40,
    color="#56B4E9",
    edgecolor="black",
)

axes[1].axvline(
    weighted_mean_depth,
    color="red",
    linestyle="--",
    label=f"Weighted mean: {weighted_mean_depth:.2f}×",
)

axes[1].set_xlabel("Mean depth per contig")
axes[1].set_ylabel("Number of contigs")
axes[1].set_title("Contig depth distribution")
axes[1].legend()

axes[1].set_xlim(0, 50)

axes[1].text(
    0.98,
    0.85,
    f"{depth_above_50_contigs} contigs above 50×",
    transform=axes[1].transAxes,
    horizontalalignment="right",
)

plt.tight_layout()
plt.savefig(figure_file, dpi=300, bbox_inches="tight")

print(summary.to_string(index=False))
print(f"\nSaved: {summary_file}")
print(f"Saved: {flagged_file}")
print(f"Saved: {figure_file}")
