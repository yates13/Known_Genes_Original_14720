#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=06.5b_mosq_cov
#SBATCH --output=../logs/06.5b_mosq_cov_%j.out
#SBATCH --error=../logs/06.5b_mosq_cov_%j.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu
# 06.5b - rerun of the mosquito DIAMOND searches with subject length and coverage columns.
# Cquinq/Aedes: all unlabeled genes (unlabeled_genes.fa). Agam/Dmel: genes still unlabeled after
# those two (still_unlabeled.fa), as in the original Oct 15 run. Writes to results/mosquito_diamond_cov/;
# the original tables are not touched. Run from scripts/.
set -eo pipefail
source /scratch/alpine/c832500103@colostate.edu/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/alpine/c832500103@colostate.edu/blast_env
set -u
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
DIAMOND=/scratch/alpine/c832500103@colostate.edu/conda_envs/annot_env/bin/diamond
[ -x "$DIAMOND" ] || { echo "diamond missing" >&2; exit 1; }
"$DIAMOND" version
OLD=../results/mosquito_diamond
OUT=../results/mosquito_diamond_cov
mkdir -p $OUT

# db  query_fasta  evalue  mode   (Agam/Dmel settings are from the Oct 15 command history)
while read -r S IN EVAL MODE; do
  [ -s "$IN" ] || { echo "missing $IN" >&2; exit 1; }
  [ -f ../results/diamond_db/$S.dmnd ] || "$DIAMOND" makedb --in ../results/mosquito_proteins/$S.faa -d ../results/diamond_db/$S < /dev/null
  "$DIAMOND" blastx -d ../results/diamond_db/$S -q "$IN" \
    --evalue "$EVAL" --max-target-seqs 5 "$MODE" --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle slen sstart send qcovhsp scovhsp \
    -o $OUT/$S.tsv.tmp < /dev/null
  mv $OUT/$S.tsv.tmp $OUT/$S.tsv
  echo "$S: queries $(grep -c '>' $IN); genes with a hit $(cut -f1 $OUT/$S.tsv | sort -u | wc -l) (old: $(cut -f1 $OLD/$S.tsv | sort -u | wc -l))"
  echo "$S: rows differing from the old table in columns 1-7: $(diff <(cut -f1-7 $OUT/$S.tsv | sort) <(sort $OLD/$S.tsv) | grep -c '^[<>]' || true)"
done <<EOF
Cquinq $OLD/unlabeled_genes.fa 1e-10 --sensitive
Aedes $OLD/unlabeled_genes.fa 1e-10 --sensitive
Agam $OLD/still_unlabeled.fa 1e-5 --ultra-sensitive
Dmel $OLD/still_unlabeled.fa 1e-5 --ultra-sensitive
EOF
