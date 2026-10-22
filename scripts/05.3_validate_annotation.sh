#!/bin/bash
# 05.3_validate_annotation.sh - sanity checks on the final v3 GTF (run from scripts/)
R=../results
V=$R/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.gtf
G=$R/basefeatures_realgenes.annotated.mosquito_rescue.gtf
echo "lines (expect 14720 / 14720):"; wc -l < $G; wc -l < $V
echo "cols 1-8 diff vs rescue GTF (expect 0):"; diff <(cut -f1-8 $G) <(cut -f1-8 $V) | wc -l
echo "unannotated (expect 1694):"; grep -c 'annotation_source "none"' $V
echo "confidence counts (expect 9546 high / 561 low / 2919 medium, no unverified):"
grep -o 'annotation_confidence "[^"]*"' $V | sort | uniq -c
echo "uncharacterized protein (expect 561):"; grep -c 'annotation_name "uncharacterized protein"' $V
echo "LOW QUALITY names left (expect 0):"; grep -c 'annotation_name "LOW QUALITY' $V
echo "species (expect 4 species, 2227 total, no unknown):"
grep -o 'annotation_species "[^"]*"' $V | sort | uniq -c
echo "gene5028 (expect eL32):"; grep 'gene_id "gene5028";' $V | grep -o 'annotation_name "[^"]*"'
echo "gene9910 (expect demoted):"; grep 'gene_id "gene9910";' $V | grep -o 'annotation_name "[^"]*"'
