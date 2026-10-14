#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --job-name=06.5_mosq_diamond
#SBATCH --output=../logs/06.5_mosq_diamond_%j.out
#SBATCH --error=../logs/06.5_mosq_diamond_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu
# 06.5 - DIAMOND blastx of the genes still unlabeled after SwissProt (annotation_source "none")
# against C. quinquefasciatus then Aedes RefSeq proteins. Run from scripts/ (sbatch or bash).
set -eo pipefail
source /scratch/alpine/c832500103@colostate.edu/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/alpine/c832500103@colostate.edu/blast_env
set -u
cd "${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
DIAMOND=/scratch/alpine/c832500103@colostate.edu/conda_envs/annot_env/bin/diamond
[ -x "$DIAMOND" ] || { echo "diamond missing" >&2; exit 1; }
SAMTOOLS=/scratch/alpine/$USER/conda_envs/hisat2_env/bin/samtools
OUT=../results/mosquito_diamond
mkdir -p $OUT

# 1. unlabeled gene IDs from the current annotated GTF
awk -F'\t' '$3=="gene" && $9 ~ /annotation_source "none"/{match($9,/gene_id "[^"]+"/); print substr($9,RSTART+9,RLENGTH-10)}' \
  ../results/basefeatures_realgenes.annotated.gtf > $OUT/unlabeled_gene_ids.txt
echo "unlabeled genes: $(wc -l < $OUT/unlabeled_gene_ids.txt)"

# 2. their per-gene sequences (same known_genes.fa used for SwissProt)
$SAMTOOLS faidx ../results/known_genes.fa -r $OUT/unlabeled_gene_ids.txt > $OUT/unlabeled_genes.fa
echo "sequences extracted: $(grep -c '>' $OUT/unlabeled_genes.fa)"

# 3. DIAMOND databases (built once) and search; e-value 1e-10, top 5 targets
for S in Cquinq Aedes; do
  [ -f ../results/diamond_db/$S.dmnd ] || "$DIAMOND" makedb --in ../results/mosquito_proteins/$S.faa -d ../results/diamond_db/$S
  "$DIAMOND" blastx -d ../results/diamond_db/$S -q $OUT/unlabeled_genes.fa \
    --evalue 1e-10 --max-target-seqs 5 --sensitive --threads "${SLURM_CPUS_PER_TASK:-4}" \
    --outfmt 6 qseqid sseqid pident length evalue bitscore stitle -o $OUT/$S.tsv.tmp
  mv $OUT/$S.tsv.tmp $OUT/$S.tsv
  echo "$S: $(cut -f1 $OUT/$S.tsv | sort -u | wc -l) genes with a hit"
done

