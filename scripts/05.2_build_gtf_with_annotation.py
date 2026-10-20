#!/usr/bin/env python3
"""05.2 - write a GTF whose gene lines record where each annotation came from.
Evidence files are read in priority order; the first file with a hit for a gene wins.
gene_id is never modified. Only feature == gene lines are changed. Run from scripts/."""
import re
import sys
from collections import Counter

IN_GTF = "../data/basefeatures_realgenes.gtf"
OUT_GTF = "../results/basefeatures_realgenes.annotated.gtf"
OUT_LOOKUP = "../results/annotation_lookup.tsv"
EVIDENCE = [  # priority order - add new evidence files here
    "../results/annotation_evidence.SwissProt_BLASTX.best.tsv"
    "../results/annotation_evidence.Mosquito_DIAMOND.filtered.tsv",
]

GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')


def clean(v):
    return v.replace('"', "'").replace(";", ",").replace("\t", " ").replace("\n", " ").strip()


evidence = {}  # gene_id -> rows, in priority order
for path in EVIDENCE:
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(cols, line.rstrip("\n").split("\t")))
            evidence.setdefault(row["gene_id"], []).append(row)


def label(r):
    return f'{r["source"]}_{r["method"]}'


def annotation_attrs(rows):
    w = rows[0]
    pairs = [
        ("annotation_source", w["source"]), ("annotation_method", w["method"]),
        ("annotation_accession", w["accession"]), ("annotation_name", w["name"]),
        ("annotation_evalue", w["evalue"]), ("annotation_bitscore", w["bitscore"]),
        ("annotation_pident", w["pident"]), ("annotation_tier", w["tier"]),
        ("annotation_transcript", w["transcript_id"]),
        ("annotation_shared_transcript", w["genes_sharing_transcript"]),
    ]
    others = []
    for r in rows[1:]:
        if label(r) != label(w) and label(r) not in others:
            others.append(label(r))
    if others:
        pairs.append(("annotation_other_sources", ",".join(others)))
    # same format as the original 04.6 note, kept for downstream compatibility
    pairs.append(("note", f'{w["name"]} [{w["method"].lower()}_{w["source"].lower()}: {w["tier"]}]'))
    return "".join(f' {k} "{clean(v)}";' for k, v in pairs)


counts = Counter()
seen_genes = set()
with open(IN_GTF) as fin, open(OUT_GTF, "w") as fout, open(OUT_LOOKUP, "w") as flk:
    flk.write("gene_id\tsource\tname\ttier\n")
    for line in fin:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9 or cols[2] != "gene":
            fout.write(line)
            continue
        m = GENE_ID_RE.search(cols[8])
        gid = m.group(1) if m else None
        attrs = cols[8].rstrip()
        if not attrs.endswith(";"):
            attrs += ";"
        rows = evidence.get(gid) if gid else None
        if rows:
            attrs += annotation_attrs(rows)
            counts[label(rows[0])] += 1
            flk.write(f'{gid}\t{label(rows[0])}\t{rows[0]["name"]}\t{rows[0]["tier"]}\n')
        else:
            attrs += ' annotation_source "none";'
            counts["none"] += 1
            if gid:
                flk.write(f"{gid}\tnone\t\t\n")
        if gid:
            seen_genes.add(gid)
        cols[8] = attrs
        fout.write("\t".join(cols) + "\n")

print(f"gene lines by annotation source: {dict(counts)}", file=sys.stderr)
print(f"evidence genes NOT found in the GTF: {len(set(evidence) - seen_genes)}", file=sys.stderr)
print(f"wrote {OUT_GTF}", file=sys.stderr)
print(f"wrote {OUT_LOOKUP}", file=sys.stderr)
