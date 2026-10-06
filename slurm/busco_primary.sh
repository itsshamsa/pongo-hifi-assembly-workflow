#!/bin/bash -l

#SBATCH --job-name=busco_primary
#SBATCH --partition=iota
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --output=logs/busco_primary_%j.out
#SBATCH --error=logs/busco_primary_%j.err


set -euo pipefail


echo "Started: $(date)"
echo "Running on: $(hostname)"


busco \
    -i assembly_fasta/pongo_denovo.primary.fa \
    -o primary_primates_odb10_busco \
    -m genome \
    -l primates_odb10 \
    -c "$SLURM_CPUS_PER_TASK" \
    --out_path results/busco \
    --download_path busco_downloads


echo "Finished: $(date)"
