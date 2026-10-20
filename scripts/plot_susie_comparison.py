import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


parser = argparse.ArgumentParser(
    description="Plot Susie assembly comparison results."
)

parser.add_argument(
    "--assembly-stats",
    required=True,
)

parser.add_argument(
    "--alignment-summary",
    required=True,
)

parser.add_argument(
    "--output",
    required=True,
)

args = parser.parse_args()


assembly_stats = pd.read_csv(
    args.assembly_stats,
    sep="\t",
)

alignment_summary = pd.read_csv(
    args.alignment_summary,
    sep="\t",
)


assembly_names = {
    "assembly_fasta/pongo_denovo.primary.fa": "Our primary",
    "assembly_fasta/pongo_denovo.hap1.fa": "Our hap1",
    "assembly_fasta/pongo_denovo.hap2.fa": "Our hap2",
    "reference_assemblies/susie_ncbi/primary.fasta":
        "NCBI primary",
    "reference_assemblies/susie_ncbi/alternate.fasta":
        "NCBI alternate",
}

assembly_stats["label"] = assembly_stats["file"].map(
    assembly_names
)

assembly_stats["n50_mb"] = (
    assembly_stats["N50"] / 1_000_000
)


comparison_names = {
    "our_hap1_vs_ncbi_primary":
        "Hap1 → NCBI primary",
    "our_hap1_vs_ncbi_alternate":
        "Hap1 → NCBI alternate",
    "our_hap2_vs_ncbi_primary":
        "Hap2 → NCBI primary",
    "our_hap2_vs_ncbi_alternate":
        "Hap2 → NCBI alternate",
    "our_primary_vs_ncbi_primary":
        "Primary → NCBI primary",
    "our_primary_vs_ncbi_alternate":
        "Primary → NCBI alternate",
}

alignment_summary["label"] = (
    alignment_summary["comparison"].map(
        comparison_names
    )
)


colors = [
    "#0072B2",
    "#56B4E9",
    "#56B4E9",
    "#D55E00",
    "#E69F00",
]

figure, axes = plt.subplots(
    2,
    2,
    figsize=(15, 11),
)


# Panel A: number of contigs
axes[0, 0].bar(
    assembly_stats["label"],
    assembly_stats["num_seqs"],
    color=colors,
)

axes[0, 0].set_yscale("log")
axes[0, 0].set_ylabel("Number of contigs (log scale)")
axes[0, 0].set_title("A. Assembly fragmentation")
axes[0, 0].tick_params(axis="x", rotation=30)


# Panel B: N50
axes[0, 1].bar(
    assembly_stats["label"],
    assembly_stats["n50_mb"],
    color=colors,
)

axes[0, 1].set_ylabel("N50 (Mb)")
axes[0, 1].set_title("B. Assembly contiguity")
axes[0, 1].tick_params(axis="x", rotation=30)


# Panel C: alignment breadth
positions = range(len(alignment_summary))
bar_height = 0.35

axes[1, 0].barh(
    [position - bar_height / 2 for position in positions],
    alignment_summary["query_coverage_percent"],
    height=bar_height,
    label="Our assembly aligned",
    color="#0072B2",
)

axes[1, 0].barh(
    [position + bar_height / 2 for position in positions],
    alignment_summary["reference_coverage_percent"],
    height=bar_height,
    label="NCBI reference covered",
    color="#E69F00",
)

axes[1, 0].set_yticks(list(positions))
axes[1, 0].set_yticklabels(
    alignment_summary["label"]
)

axes[1, 0].set_xlim(90, 100)
axes[1, 0].set_xlabel("Aligned sequence (%)")
axes[1, 0].set_title("C. Whole-genome alignment breadth")
axes[1, 0].legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.18),
    ncol=2,
    frameon=False,
)
axes[1, 0].invert_yaxis()


# Panel D: estimated identity
axes[1, 1].scatter(
    alignment_summary[
        "estimated_identity_percent"
    ],
    list(positions),
    color="#009E73",
    s=80,
)

axes[1, 1].set_yticks(list(positions))
axes[1, 1].set_yticklabels(
    alignment_summary["label"]
)

axes[1, 1].invert_yaxis()

axes[1, 1].set_xlabel(
    "Estimated nucleotide identity (%)"
)

axes[1, 1].set_title(
    "D. Sequence agreement"
)

identity_min = (
    alignment_summary[
        "estimated_identity_percent"
    ].min()
)

identity_max = (
    alignment_summary[
        "estimated_identity_percent"
    ].max()
)

axes[1, 1].set_xlim(
    identity_min - 0.005,
    identity_max + 0.005,
)

axes[1, 1].grid(
    axis="x",
    linestyle="--",
    alpha=0.4,
)


figure.suptitle(
    "Susie assembly comparison: "
    "new de novo assembly versus NCBI",
    fontsize=16,
)

figure.tight_layout()

output_path = Path(args.output)
output_path.parent.mkdir(
    parents=True,
    exist_ok=True,
)

figure.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
)

print(f"Saved: {output_path}")
