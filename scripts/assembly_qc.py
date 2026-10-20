import argparse
from pathlib import Path

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


def main():
    parser = argparse.ArgumentParser(
        description="Calculate assembly statistics from three FASTA files."
    )

    parser.add_argument(
        "--primary",
        required=True,
        help="Path to the primary assembly FASTA."
    )

    parser.add_argument(
        "--hap1",
        required=True,
        help="Path to the haplotype 1 FASTA."
    )

    parser.add_argument(
        "--hap2",
        required=True,
        help="Path to the haplotype 2 FASTA."
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Path for the output TSV summary."
    )

    args = parser.parse_args()

    assemblies = {
        "primary": args.primary,
        "hap1": args.hap1,
        "hap2": args.hap2
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

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    print(results_table.to_string(index=False))

    results_table.to_csv(
        output_path,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()
