#!/usr/bin/env python3
"""
06.7_build_mosquito_rescue_gtf.py

Adds filtered mosquito DIAMOND annotations to genes that remain
unannotated in the SwissProt-derived GTF.

Input:
    ../results/basefeatures_realgenes.annotated.gtf
#    ../results/annotation_evidence.Mosquito_DIAMOND.filtered.tsv
    ../results/annotation_evidence.Mosquito_DIAMOND.tsv

Output:
    ../results/basefeatures_realgenes.annotated.mosquito_rescue.gtf

Only genes with:

    annotation_source "none"

are modified.

All existing SwissProt annotations are preserved.
"""

import csv
import re

INPUT_GTF = "../results/basefeatures_realgenes.annotated.gtf"

#old Filter with quin and aed only
#INPUT_EVIDENCE = "../results/annotation_evidence.Mosquito_DIAMOND.filtered.tsv"

INPUT_EVIDENCE = "../results/annotation_evidence.Mosquito_DIAMOND.tsv"
OUTPUT_GTF = "../results/basefeatures_realgenes.annotated.mosquito_rescue.gtf"

print("Loading mosquito evidence...")

mosquito_hits = {}

with open(INPUT_EVIDENCE) as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for row in reader:
        mosquito_hits[row["gene_id"]] = row

print(f"Loaded {len(mosquito_hits):,} mosquito annotations")

updated = 0
unchanged = 0

with open(INPUT_GTF) as fin, open(OUTPUT_GTF, "w") as fout:

    for line in fin:

        if line.startswith("#"):
            fout.write(line)
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 9:
            fout.write(line)
            continue

        if fields[2] != "gene":
            fout.write(line)
            continue

        attrs = fields[8]

        m = re.search(r'gene_id "([^"]+)"', attrs)

        if not m:
            fout.write(line)
            continue

        gene_id = m.group(1)

        if gene_id not in mosquito_hits:
            fout.write(line)
            continue

        if 'annotation_source "none"' not in attrs:
            fout.write(line)
            continue

        hit = mosquito_hits[gene_id]
        safe_name = hit["name"].replace('"', "'").replace(";", ",")
        replacement = (
            f'annotation_source "{hit["source"]}"; '
            f'annotation_method "{hit["method"]}"; '
            f'annotation_accession "{hit["accession"]}"; '
            f'annotation_name "{safe_name}"; '
            f'annotation_evalue "{hit["evalue"]}"; '
            f'annotation_bitscore "{hit["bitscore"]}"; '
            f'annotation_pident "{hit["pident"]}"; '
            f'annotation_tier "{hit["tier"]}"; '
            f'annotation_transcript "{hit["transcript_id"]}"; '
            f'annotation_shared_transcript "{hit["genes_sharing_transcript"]}"; '
            f'note "{safe_name} [diamond_mosquito_refseq: {hit["tier"]}]";'
        )
        attrs = attrs.replace('annotation_source "none";', replacement)
        fields[8] = attrs
        fout.write("\t".join(fields) + "\n")
        updated += 1

print()
print(f"Genes updated:   {updated:,}")
print(f"Output written:  {OUTPUT_GTF}")
