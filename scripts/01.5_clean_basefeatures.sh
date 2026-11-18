#!/bin/bash
#SBATCH --partition=acpu
#SBATCH --job-name=01.5_clean_basefeatures
#SBATCH --output=../logs/01.5_clean_basefeatures_%j.out
#SBATCH --error=../logs/01.5_clean_basefeatures_%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --qos=cpu-normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=c832500103@colostate.edu


# STEP 01.5 (from Culex_Tarsalis_GTF 01_novel_locus_pipeline/01.5_Real_Contigs.sh) - strip the raw basefeatures GTF down to real gene records on
# real contigs only, with Mt contig naming fixed BEFORE filtering so the
# mitochondrial genes survive. Independent of the GMAP steps.
# Submit from inside scripts/.

set -euo pipefail

# 0. fix Mt contig naming BEFORE any filtering - the raw GTF uses plain "Mt"
#    but the genome FASTA uses "Mt_CtarK1_noOVERLAP_RC.1". If this isn't done
#    first, the contig filter below silently drops all 6 mitochondrial genes
#    since "Mt" never matches anything in real_contigs.txt.
sed 's/^Mt\t/Mt_CtarK1_noOVERLAP_RC.1\t/' ../data/basefeatures_raw.gtf > ../results/basefeatures_raw.mtfixed.gtf
echo "Mt rows after rename (expect 47 - gene/transcript/exon/CDS/UTR lines for 6 genes):"
grep -c "Mt_CtarK1_noOVERLAP_RC.1" ../results/basefeatures_raw.mtfixed.gtf

# 1. list of contigs actually present in the genome build
cut -f1 ../data/genome.fa.fai > ../results/real_contigs.txt   # contig names from the genome index (genome.fa was purged from scratch)
echo "Real contigs found:"
wc -l ../results/real_contigs.txt

# 2. drop rows on any contig not in real_contigs.txt
awk -F'\t' 'NR==FNR{keep[$1]=1; next} ($1 in keep)' \
    ../results/real_contigs.txt ../results/basefeatures_raw.mtfixed.gtf > ../results/basefeatures_clean.gtf
echo "After contig filter:"
wc -l ../results/basefeatures_clean.gtf

# 3. keep only real gene predictions - source column (col2) must be maker,
#    Geneious, or VectorBase. Everything else tagged "gene" in col3
#    (repeatmasker, protein2genome, blastx, est2genome, blastn, repeatrunner)
#    is an evidence/repeat track, not a real gene call.
awk -F'\t' '$3=="gene" && ($2=="maker" || $2=="Geneious" || $2=="VectorBase")' \
    ../results/basefeatures_clean.gtf > ../results/basefeatures_realgenes.gtf
echo "Real gene records (expect 14720 exactly: 14713 maker + 6 Geneious + 1 VectorBase):"
wc -l ../results/basefeatures_realgenes.gtf
echo "Mt genes present (expect 6):"
grep -c "Mt_CtarK1_noOVERLAP_RC.1" ../results/basefeatures_realgenes.gtf
