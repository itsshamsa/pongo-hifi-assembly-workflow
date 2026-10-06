#!/bin/bash -l

#SBATCH --job-name=busco_hap2
#SBATCH --partition=iota
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --output=logs/busco_hap2_%j.out
#SBATCH --error=logs/busco_hap2_%j.err

set -euo pipefail

echo "Started: $(date)"
echo "Running on: $(hostname)"

busco \
    -i assembly_fasta/pongo_denovo.hap2.fa \
    -o hap2_primates_odb10_busco \
    -m genome \
    -l primates_odb10 \
    -c "$SLURM_CPUS_PER_TASK" \
    --out_path results/busco \
    --download_path busco_downloads

echo "Finished: $(date)"
