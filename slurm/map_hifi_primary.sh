#!/bin/bash -l

#SBATCH --job-name=map_primary
#SBATCH --partition=iota
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=logs/map_primary_%j.out
#SBATCH --error=logs/map_primary_%j.err

set -euo pipefail

REFERENCE="assembly_fasta/pongo_denovo.primary.fa"
READS="full.fastq.gz"
BAM="alignments/pongo_denovo.primary.hifi.sorted.bam"
TEMP_PREFIX="alignments/tmp/primary_sort"

echo "Started: $(date)"
echo "Running on: $(hostname)"
echo "Reference: $REFERENCE"
echo "Reads: $READS"

minimap2 \
    -ax map-hifi \
    -t 10 \
    "$REFERENCE" \
    "$READS" |
samtools sort \
    -@ 5 \
    -m 2G \
    -T "$TEMP_PREFIX" \
    -o "$BAM" \
    -

samtools index -@ 10 "$BAM"

samtools flagstat -@ 10 "$BAM" \
    > results/mapping/primary.flagstat.txt

echo "BAM: $BAM"
echo "Finished: $(date)"
