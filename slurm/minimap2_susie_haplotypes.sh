#!/bin/bash -l

#SBATCH --job-name=susie_asm
#SBATCH --partition=iota
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=12:00:00
#SBATCH --array=0-3%2
#SBATCH --output=logs/minimap2_susie_%A_%a.out
#SBATCH --error=logs/minimap2_susie_%A_%a.err

set -euo pipefail

NCBI_DIR="reference_assemblies/susie_ncbi"
OUR_DIR="assembly_fasta"
RESULT_DIR="results/susie_comparison/minimap2"

mkdir -p "$RESULT_DIR"

case "$SLURM_ARRAY_TASK_ID" in
    0)
        REFERENCE="$NCBI_DIR/primary.fasta"
        QUERY="$OUR_DIR/pongo_denovo.hap1.fa"
        NAME="our_hap1_vs_ncbi_primary"
        ;;
    1)
        REFERENCE="$NCBI_DIR/alternate.fasta"
        QUERY="$OUR_DIR/pongo_denovo.hap1.fa"
        NAME="our_hap1_vs_ncbi_alternate"
        ;;
    2)
        REFERENCE="$NCBI_DIR/primary.fasta"
        QUERY="$OUR_DIR/pongo_denovo.hap2.fa"
        NAME="our_hap2_vs_ncbi_primary"
        ;;
    3)
        REFERENCE="$NCBI_DIR/alternate.fasta"
        QUERY="$OUR_DIR/pongo_denovo.hap2.fa"
        NAME="our_hap2_vs_ncbi_alternate"
        ;;

    4)
        REFERENCE="$NCBI_DIR/primary.fasta"
        QUERY="$OUR_DIR/pongo_denovo.primary.fa"
        NAME="our_primary_vs_ncbi_primary"
        ;;
    5)
        REFERENCE="$NCBI_DIR/alternate.fasta"
        QUERY="$OUR_DIR/pongo_denovo.primary.fa"
        NAME="our_primary_vs_ncbi_alternate"
        ;;

esac

echo "Started: $(date)"
echo "Reference: $REFERENCE"
echo "Query: $QUERY"
echo "Comparison: $NAME"

minimap2 \
    -x asm5 \
    --secondary=no \
    -t "$SLURM_CPUS_PER_TASK" \
    "$REFERENCE" \
    "$QUERY" \
    > "$RESULT_DIR/$NAME.paf"

echo "Finished: $(date)"
