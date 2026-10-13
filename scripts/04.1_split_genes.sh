#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p ../data/chunks ../logs
rm -f ../data/chunks/known_genes.*.fa
awk -v n=20 '/^>/{i=(c++)%n} {print > sprintf("../data/chunks/known_genes.%02d.fa", i+1)}' ../results/known_genes.fa
for f in ../data/chunks/known_genes.*.fa; do echo "$f $(grep -c '>' $f)"; done
echo "total_records $(cat ../data/chunks/known_genes.*.fa | grep -c '>')"
echo "duplicate_ids $(grep -h '>' ../data/chunks/known_genes.*.fa | sort | uniq -d | wc -l)"

