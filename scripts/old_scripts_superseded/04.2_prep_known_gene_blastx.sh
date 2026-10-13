#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=04.2_prep_known_gene_blastx
#SBATCH --output=../logs/04.2_prep_known_gene_blastx_%j.out
#SBATCH --error=../logs/04.2_prep_known_gene_blastx_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --qos=cpu-normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu
#
# STEP 04.2 - extract representative-transcript sequences for the 14,720
# known genes (13,953 with a supporting transcript, from
# gene_to_representative_transcript.tsv / representative_transcript_ids_knowngenes.txt),
# then chunk for a BLASTX array job. Mirrors 01.13 steps 5-6 exactly, so the
# known-gene BLASTX run is directly comparable to the novel-locus one.
#
# Depends on: 01_build_known_gene_representatives.py (already run)
set -euo pipefail

# 1. pull sequences out of the transcriptome (same tool/env as 01.13)
source /scratch/alpine/c832500103@colostate.edu/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/cellSquito

seqkit grep -f ../results/representative_transcript_ids_knowngenes.txt ../data/transcriptome.fasta \
    > ../results/known_gene_representative.fasta

echo "expected sequences:"
wc -l ../results/representative_transcript_ids_knowngenes.txt
echo "sequences actually pulled:"
grep -c "^>" ../results/known_gene_representative.fasta

# 2. split into chunks for the BLASTX array job
#    13,953 sequences vs 11,993 for the novel set - same order of magnitude,
#    so reuse 20 chunks for a similar per-chunk runtime.
mkdir -p ../results/known_gene_chunks
seqkit split2 -p 20 -O ../results/known_gene_chunks ../results/known_gene_representative.fasta

echo "chunk files created:"
ls ../results/known_gene_chunks/
