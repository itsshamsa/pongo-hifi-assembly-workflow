#!/bin/bash -l

#SBATCH --job-name=seqkit_stats
#SBATCH --partition=iota
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=logs/seqkit_%j.out
#SBATCH --error=logs/seqkit_%j.err


set -euo pipefail


echo "Started: $(date)"
echo "Running on: $(hostname)"


seqkit stats \
    -a \
    -T \
    -j "$SLURM_CPUS_PER_TASK" \
    assembly_fasta/*.fa \
    > results/seqkit_stats.tsv


echo "Finished: $(date)"
