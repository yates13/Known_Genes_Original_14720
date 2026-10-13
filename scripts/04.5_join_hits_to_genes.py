#!/usr/bin/env python3
"""04.5 - join per-GENE BLASTX best hits to the known genes.
Queries were each gene's own sequence (known_genes.fa), so qseqid IS the gene_id.
Output columns match the old version (transcript_id = gene_id, genes_sharing_transcript = 1)
so 05.1 and 05.2 run unchanged. Run from scripts/."""
import sys
from collections import Counter

FA = "../results/known_genes.fa"
HITS = "../results/known_gene_blastx_best_hits_parsed.tsv"
OUT = "../results/known_gene_annotation_by_gene.tsv"

hits = {}
with open(HITS) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6:
            continue
        try:
            float(f[1])          # skips the header row
        except ValueError:
            continue
        hits[f[0]] = f[1:6]

genes = []
with open(FA) as fh:
    for line in fh:
        if line.startswith(">"):
            genes.append(line[1:].split()[0])
gene_set = set(genes)

tiers = Counter()
n_hit = 0
with open(OUT, "w") as out:
    out.write("gene_id\ttranscript_id\tgenes_sharing_transcript\thas_hit\tevalue\tbitscore\ttier\tname\tstitle\n")
    for g in genes:
        h = hits.get(g)
        if h:
            n_hit += 1
            tiers[h[2]] += 1
            out.write(f"{g}\t{g}\t1\tyes\t" + "\t".join(h) + "\n")
        else:
            out.write(f"{g}\t{g}\t1\tno\t\t\t\t\t\n")

orphans = [q for q in hits if q not in gene_set]
print(f"genes in known_genes.fa:   {len(genes)} ({len(gene_set)} unique)", file=sys.stderr)
print(f"genes with a BLASTX hit:   {n_hit} ({100*n_hit/len(genes):.1f}%)", file=sys.stderr)
print(f"hit IDs not in the fasta:  {len(orphans)}  (want 0)", file=sys.stderr)
print(f"tier breakdown:            {dict(tiers)}", file=sys.stderr)
print(f"wrote {OUT}", file=sys.stderr)
