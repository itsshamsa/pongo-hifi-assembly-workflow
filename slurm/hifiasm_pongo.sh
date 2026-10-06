#!/bin/bash -l

#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --partition=iota
#SBATCH --job-name=pongo_hifiasm
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err

echo "Started: $(date)"
echo "Running on: $(hostname)"

hifiasm -o pongo_denovo -t 10 full.fastq.gz

echo "Finished: $(date)"
