import numpy as np
import pandas as pd


def read_fasta_lengths(fasta_file):
    lengths = []
    current_length = 0

    with open(fasta_file, "r") as fasta:
        for line in fasta:
            line = line.strip()

            if line.startswith(">"):
                if current_length > 0:
                    lengths.append(current_length)

                current_length = 0

            else:
                current_length = current_length + len(line)

    if current_length > 0:
        lengths.append(current_length)

    return lengths


def calculate_n50(lengths):
    sorted_lengths = sorted(lengths, reverse=True)
    half_assembly = sum(lengths) / 2
    cumulative_length = 0

    for position, length in enumerate(sorted_lengths, start=1):
        cumulative_length = cumulative_length + length

        if cumulative_length >= half_assembly:
            return length, position


assemblies = {
    "hap1": "assembly_fasta/pongo_denovo.hap1.fa",
    "hap2": "assembly_fasta/pongo_denovo.hap2.fa",
    "primary": "assembly_fasta/pongo_denovo.primary.fa"
}

results = []


for assembly_name, fasta_file in assemblies.items():
    contig_lengths = read_fasta_lengths(fasta_file)

    total_length = sum(contig_lengths)
    n50, l50 = calculate_n50(contig_lengths)

    assembly_stats = {
        "assembly": assembly_name,
        "contigs": len(contig_lengths),
        "total_bp": total_length,
        "longest_bp": max(contig_lengths),
        "mean_bp": np.mean(contig_lengths),
        "median_bp": np.median(contig_lengths),
        "n50_bp": n50,
        "l50": l50
    }

    results.append(assembly_stats)


results_table = pd.DataFrame(results)

print(results_table.to_string(index=False))

results_table.to_csv(
    "results/assembly_summary.tsv",
    sep="\t",
    index=False
)