#!/bin/bash -l

#SBATCH --job-name=coverage_primary
#SBATCH --partition=iota
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=logs/coverage_primary_%j.out
#SBATCH --error=logs/coverage_primary_%j.err

set -euo pipefail

BAM="alignments/pongo_denovo.primary.hifi.sorted.bam"

echo "Started: $(date)"
echo "Running on: $(hostname)"

samtools idxstats "$BAM" \
    > results/mapping/primary.idxstats.tsv

samtools coverage "$BAM" \
    -o results/mapping/primary.coverage.tsv

echo "Finished: $(date)"
