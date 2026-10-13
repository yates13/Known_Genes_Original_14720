#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=04.3_known_gene_blastx_array
#SBATCH --output=../logs/04.3_known_gene_blastx_array_%A_%a.out
#SBATCH --error=../logs/04.3_known_gene_blastx_array_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --qos=cpu-normal
#SBATCH --array=1-20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu
#
# STEP 04.3 - BLASTX the known-gene representative transcripts against
# swissprot, one chunk per array task. Identical settings to 01.14 (same
# db, same outfmt, same evalue/max_target_seqs) so results are directly
# comparable/mergeable with the novel-locus blastx output.
#
# Depends on: 04.2_prep_known_gene_blastx.sh
set -euo pipefail
source /scratch/alpine/c832500103@colostate.edu/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/alpine/c832500103@colostate.edu/blast_env

mkdir -p ../results/known_gene_blastx_results

PART=$(printf "%03d" $SLURM_ARRAY_TASK_ID)
QUERY="../results/known_gene_chunks/known_gene_representative.part_${PART}.fasta"
OUT="../results/known_gene_blastx_results/blastx_part_${PART}.tsv"

if [ ! -s "$QUERY" ]; then
    echo "ERROR: query file $QUERY is missing or empty. Check 02.5 completed correctly." >&2
    exit 1
fi

blastx -query "$QUERY" \
       -db ../results/blastx_db/swissprot \
       -out "$OUT" \
       -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
       -evalue 1e-5 -max_target_seqs 1 -num_threads 8
