#!/bin/bash
#SBATCH --job-name=kg_blastx
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --array=1-20
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu
#SBATCH --output=../logs/kg_blastx_%A_%a.out
#SBATCH --error=../logs/kg_blastx_%A_%a.err
set -euo pipefail
source /scratch/alpine/c832500103@colostate.edu/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/alpine/c832500103@colostate.edu/blast_env
cd "$SLURM_SUBMIT_DIR"
IDX=$(printf "%02d" "$SLURM_ARRAY_TASK_ID")
Q=../data/chunks/known_genes.${IDX}.fa
OUT=../results/blastx_chunks/known_genes.${IDX}.blastx.tsv
mkdir -p ../results/blastx_chunks
[ -f "$Q" ] || { echo "missing $Q"; exit 1; }
blastx -query "$Q" -db ../results/blastx_db/swissprot -evalue 1e-5 -max_target_seqs 5 -max_intron_length 20000 -num_threads "$SLURM_CPUS_PER_TASK" -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" > "${OUT}.tmp"
mv "${OUT}.tmp" "$OUT"
echo "done $IDX"
