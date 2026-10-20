from collections import defaultdict
from pathlib import Path

import pandas as pd


STATS_FILE = Path(
    "results/susie_comparison/assembly_stats.tsv"
)

PAF_DIR = Path(
    "results/susie_comparison/minimap2"
)

OUTPUT_FILE = Path(
    "results/susie_comparison/minimap2_summary.tsv"
)

MIN_ALIGNMENT_LENGTH = 10_000
MIN_MAPQ = 5


comparisons = {
    "our_hap1_vs_ncbi_primary": {
        "our_haplotype": "hap1",
        "ncbi_haplotype": "primary",
        "query_file": "assembly_fasta/pongo_denovo.hap1.fa",
        "reference_file": (
            "reference_assemblies/susie_ncbi/primary.fasta"
        ),
    },
    "our_hap1_vs_ncbi_alternate": {
        "our_haplotype": "hap1",
        "ncbi_haplotype": "alternate",
        "query_file": "assembly_fasta/pongo_denovo.hap1.fa",
        "reference_file": (
            "reference_assemblies/susie_ncbi/alternate.fasta"
        ),
    },
    "our_hap2_vs_ncbi_primary": {
        "our_haplotype": "hap2",
        "ncbi_haplotype": "primary",
        "query_file": "assembly_fasta/pongo_denovo.hap2.fa",
        "reference_file": (
            "reference_assemblies/susie_ncbi/primary.fasta"
        ),
    },
    "our_hap2_vs_ncbi_alternate": {
        "our_haplotype": "hap2",
        "ncbi_haplotype": "alternate",
        "query_file": "assembly_fasta/pongo_denovo.hap2.fa",
        "reference_file": (
            "reference_assemblies/susie_ncbi/alternate.fasta"
        ),
    },

    "our_primary_vs_ncbi_primary": {
        "our_haplotype": "primary",
        "ncbi_haplotype": "primary",
        "query_file": "assembly_fasta/pongo_denovo.primary.fa",
        "reference_file": (
            "reference_assemblies/susie_ncbi/primary.fasta"
        ),
    },
    "our_primary_vs_ncbi_alternate": {
        "our_haplotype": "primary",
        "ncbi_haplotype": "alternate",
        "query_file": "assembly_fasta/pongo_denovo.primary.fa",
        "reference_file": (
            "reference_assemblies/susie_ncbi/alternate.fasta"
        ),
    },

}


def merged_interval_length(intervals_by_sequence):
    total = 0

    for intervals in intervals_by_sequence.values():
        intervals.sort()

        current_start, current_end = intervals[0]

        for start, end in intervals[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                total += current_end - current_start
                current_start, current_end = start, end

        total += current_end - current_start

    return total


def summarize_paf(paf_file):
    query_intervals = defaultdict(list)
    reference_intervals = defaultdict(list)

    alignment_count = 0
    divergence_weighted_sum = 0
    divergence_weight = 0

    with open(paf_file, "r") as paf:
        for line in paf:
            fields = line.rstrip().split("\t")

            query_name = fields[0]
            query_start = int(fields[2])
            query_end = int(fields[3])

            reference_name = fields[5]
            reference_start = int(fields[7])
            reference_end = int(fields[8])

            matches = int(fields[9])
            block_length = int(fields[10])
            mapq = int(fields[11])

            divergence = None

            for tag in fields[12:]:
                if tag.startswith("dv:f:"):
                    divergence = float(tag[5:])
                    break

            if block_length < MIN_ALIGNMENT_LENGTH:
                continue

            if mapq < MIN_MAPQ:
                continue

            alignment_count += 1

            query_intervals[query_name].append(
                (query_start, query_end)
            )

            reference_intervals[reference_name].append(
                (reference_start, reference_end)
            )

            if divergence is not None:
                divergence_weighted_sum += (
                    divergence * matches
                )
                divergence_weight += matches

    aligned_query_bp = merged_interval_length(
        query_intervals
    )

    aligned_reference_bp = merged_interval_length(
        reference_intervals
    )

    weighted_divergence = (
        divergence_weighted_sum / divergence_weight
    )

    estimated_identity = (
        1 - weighted_divergence
    ) * 100

    return {
        "alignment_blocks": alignment_count,
        "aligned_query_bp": aligned_query_bp,
        "aligned_reference_bp": aligned_reference_bp,
        "estimated_identity_percent": estimated_identity,
    }

assembly_stats = pd.read_csv(STATS_FILE, sep="\t")
assembly_lengths = dict(
    zip(assembly_stats["file"], assembly_stats["sum_len"])
)

results = []

for comparison_name, information in comparisons.items():
    paf_file = PAF_DIR / f"{comparison_name}.paf"
    summary = summarize_paf(paf_file)

    query_total = assembly_lengths[
        information["query_file"]
    ]

    reference_total = assembly_lengths[
        information["reference_file"]
    ]

    summary["comparison"] = comparison_name
    summary["our_haplotype"] = information["our_haplotype"]
    summary["ncbi_haplotype"] = information["ncbi_haplotype"]

    summary["query_coverage_percent"] = (
        summary["aligned_query_bp"] / query_total * 100
    )

    summary["reference_coverage_percent"] = (
        summary["aligned_reference_bp"]
        / reference_total
        * 100
    )

    results.append(summary)


results_table = pd.DataFrame(results)

column_order = [
    "comparison",
    "our_haplotype",
    "ncbi_haplotype",
    "alignment_blocks",
    "aligned_query_bp",
    "query_coverage_percent",
    "aligned_reference_bp",
    "reference_coverage_percent",
    "estimated_identity_percent",
]

results_table = results_table[column_order]

numeric_columns = [
    "query_coverage_percent",
    "reference_coverage_percent",
    "estimated_identity_percent",
]

results_table[numeric_columns] = (
    results_table[numeric_columns].round(4)
)

results_table.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False,
)

print(results_table.to_string(index=False))
print(f"\nSaved: {OUTPUT_FILE}")
