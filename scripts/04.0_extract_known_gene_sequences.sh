
#!/bin/bash

# 04.0 - per-gene sequences (strand-aware) for the known genes -> results/known_genes.fa

# Input : Spades_Pipeline_Organized/results/basefeatures_realgenes.gtf (14,720 gene lines)

#         Spades_Pipeline_Organized/data/genome.fa

# Output: results/known_genes.bed, results/known_genes.fa (headers = gene_id, no "(+)"/"(-)")

set -euo pipefail

cd "$(dirname "$0")"

P=/scratch/alpine/$USER/culex_tarsalis_assembly/Spades_Pipeline_Organized

ENV=/scratch/alpine/$USER/conda_envs



# 1. gene lines -> BED6 (0-based start, gene_id as name, strand kept)

awk -F'\t' '$3=="gene"{match($9,/gene_id "[^"]+"/); id=substr($9,RSTART+9,RLENGTH-10); print $1"\t"$4-1"\t"$5"\t"id"\t0\t"$7}' \

  $P/results/basefeatures_realgenes.gtf > ../results/known_genes.bed

echo "BED lines: $(wc -l < ../results/known_genes.bed)  duplicate ids: $(cut -f4 ../results/known_genes.bed | sort | uniq -d | wc -l)"



# 2. every BED contig must exist in the genome index

[ -f $P/data/genome.fa.fai ] || $ENV/hisat2_env/bin/samtools faidx $P/data/genome.fa

echo "BED contigs missing from genome (want none):"

cut -f1 ../results/known_genes.bed | sort -u | grep -vFxf <(cut -f1 $P/data/genome.fa.fai) | head || true



# 3. strand-aware sequences (LD_LIBRARY_PATH fixes the bedtools library error seen on Alpine)

LD_LIBRARY_PATH=$ENV/hisat2_env/lib $ENV/cellSquito/bin/bedtools getfasta -fi $P/data/genome.fa \

  -bed ../results/known_genes.bed -s -nameOnly -fo ../results/known_genes.fa



# 4. strip the "(+)"/"(-)" that -nameOnly -s appends, then sanity-check

sed -i 's/^\(>[^(]*\)([+-])$/\1/' ../results/known_genes.fa

echo "records: $(grep -c '>' ../results/known_genes.fa)  headers still containing '(': $(grep '>' ../results/known_genes.fa | grep -c '(' || true)  duplicate ids: $(grep '>' ../results/known_genes.fa | sort | uniq -d | wc -l)"

