#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=06.4_mosquito_blastx_pilot
#SBATCH --output=../logs/06.4_mosquito_blastx_pilot_%j.out
#SBATCH --error=../logs/06.4_mosquito_blastx_pilot_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --qos=cpu-normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu

# 06.4 - BLASTX the 5,291 SwissProt-no-hit representatives against the Aedes and Culex quinquefasciatus proteomes
set -euo pipefail
source /scratch/alpine/c832500103@colostate.edu/miniconda3/etc/profile.d/conda.sh
conda activate /scratch/alpine/c832500103@colostate.edu/blast_env

Q=../results/nohit_representatives.fasta
OUT=../results/mosquito_pilot
mkdir -p $OUT

for S in Aedes Cquinq; do
  blastx -query $Q -db ../results/mosquito_proteins/$S -evalue 1e-5 -max_target_seqs 5 -num_threads $SLURM_CPUS_PER_TASK -outfmt "6 qseqid sseqid pident length evalue bitscore stitle qlen slen qcovs" -out $OUT/$S.tsv
  cut -f1 $OUT/$S.tsv | sort -u > $OUT/$S.queries
  sort -t$'\t' -k1,1 -k5,5g $OUT/$S.tsv | awk -F'\t' '!s[$1]++' > $OUT/$S.best.tsv
  echo "$S: queries with a hit: $(wc -l < $OUT/$S.queries) of $(grep -c '^>' $Q)"
  echo "$S: best hits that are uncharacterized/hypothetical: $(grep -c -i -E 'uncharacterized|hypothetical' $OUT/$S.best.tsv)"
done

cat $OUT/Aedes.queries $OUT/Cquinq.queries | sort -u > $OUT/either.queries
echo "either proteome: $(wc -l < $OUT/either.queries)"
echo "both proteomes: $(comm -12 $OUT/Aedes.queries $OUT/Cquinq.queries | wc -l)"

