#!/usr/bin/env python3
"""05.1 - gene-level SwissProt/BLASTX results -> evidence table (one row per gene per source).
Every field comes from the single HSP that 04.5 selected (matched on evalue+bitscore), so the row is
internally consistent. Adds accession, percent identity and alignment length. Run from scripts/."""
import math
import sys
GENE_TABLE = "../results/known_gene_annotation_by_gene.tsv"
ALL_HITS = "../results/known_gene_blastx_all_hits.tsv"
OUT = "../results/annotation_evidence.SwissProt_BLASTX.tsv"
SOURCE, METHOD = "SwissProt", "BLASTX"
COLS = ["gene_id", "source", "method", "accession", "name", "evalue", "bitscore",
        "pident", "aln_len", "tier", "transcript_id", "genes_sharing_transcript"]
def accession(sseqid):
    p = sseqid.split("|")
    return p[1] if len(p) >= 3 else sseqid
# all HSPs per query (cols: qseqid sseqid pident length evalue bitscore stitle)
hsps = {}
with open(ALL_HITS) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 7:
            hsps.setdefault(f[0], []).append(f)
multi_subject = sum(1 for v in hsps.values() if len({h[1] for h in v}) > 1)
n = unmatched = 0
with open(GENE_TABLE) as fh, open(OUT, "w") as out:
    next(fh)  # header
    out.write("\t".join(COLS) + "\n")
    for line in fh:
        f = line.rstrip("\n").split("\t") + [""] * 9
        gene, tx, shared, has_hit, ev, bits, tier, name = f[:8]
        if has_hit != "yes":
            continue
        pick = None
        for h in hsps.get(tx, []):
            if math.isclose(float(h[4]), float(ev), rel_tol=1e-6) and \
               math.isclose(float(h[5]), float(bits), rel_tol=1e-6):
                pick = h
                break
        if pick is None:
            unmatched += 1
            continue
        out.write("\t".join([gene, SOURCE, METHOD, accession(pick[1]), name, ev, bits,
                             pick[2], pick[3], tier, tx, "yes" if int(shared) > 1 else "no"]) + "\n")
        n += 1
print(f"evidence rows written: {n}", file=sys.stderr)
print(f"genes with no matching HSP (not written): {unmatched}", file=sys.stderr)
print(f"queries with more than one distinct subject: {multi_subject}", file=sys.stderr)
print(f"wrote {OUT}", file=sys.stderr)
