#!/usr/bin/env python3
"""
STEP 04.5 - collapse each blastx hit file down to one best hit per query
transcript, parse a clean name out of the SwissProt stitle, and assign a
confidence tier. Run once for the known-gene results and once for the
novel-gene results (same logic, different input/output paths - see main()).

Input columns (blastx -outfmt "6 qseqid sseqid pident length evalue bitscore stitle"):
    0 qseqid   1 sseqid   2 pident   3 length   4 evalue   5 bitscore   6 stitle

"Best hit per query" = lowest evalue, tie-broken by highest bitscore.
(A query can appear multiple times because -max_target_seqs 1 still reports
multiple HSPs/segments for the same single best subject hit.)

Name parsing: SwissProt stitle looks like
    RecName: Full=Probable RNA-directed DNA polymerase ...; AltName: Full=... [Drosophila melanogaster]
We take the first "RecName: Full=...;" segment as the clean name. If no
RecName is present, fall back to the first "Full=...;" segment. If neither
semicolon-delimited form is found, we fall back to the whole stitle with any
trailing "[Organism]" tag stripped, so nothing is silently dropped.

Tiers mirror the ones already used for the novel-locus BLASTX
(strong/very strong/extremely strong), applied here to BOTH known and novel
results so they're on the same scale:
    e < 1e-40  -> extremely_strong
    e < 1e-20  -> very_strong
    e < 1e-10  -> strong
    e < 1e-5   -> any_hit
    (anything looser shouldn't appear - search was already run at -evalue 1e-5)

Output: one TSV with columns
    qseqid  evalue  bitscore  tier  parsed_name  raw_stitle
"""

import re
import sys

RECNAME_RE = re.compile(r"RecName:\s*Full=([^;]+);?")
FULL_FALLBACK_RE = re.compile(r"Full=([^;]+);?")
ORGANISM_TAIL_RE = re.compile(r"\s*\[[^\[\]]+\]\s*$")


def parse_name(stitle: str) -> str:
    m = RECNAME_RE.search(stitle)
    if m:
        return ORGANISM_TAIL_RE.sub("", m.group(1)).strip()
    m = FULL_FALLBACK_RE.search(stitle)
    if m:
        return ORGANISM_TAIL_RE.sub("", m.group(1)).strip()
    # last resort: whole stitle minus a trailing "[Organism]" tag
    return ORGANISM_TAIL_RE.sub("", stitle).strip()


def tier_for_evalue(evalue: float) -> str:
    if evalue < 1e-40:
        return "extremely_strong"
    if evalue < 1e-20:
        return "very_strong"
    if evalue < 1e-10:
        return "strong"
    if evalue < 1e-5:
        return "any_hit"
    return "below_threshold"  # shouldn't occur, search cutoff was already 1e-5


def process(in_path: str, out_path: str, label: str):
    # best[qseqid] = (evalue, -bitscore, sseqid, stitle)   (lower tuple = better; negate bitscore so higher bitscore sorts first on tie)
    best = {}
    total_lines = 0
    skipped = 0

    with open(in_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            total_lines += 1
            cols = line.split("\t")
            if len(cols) < 7:
                skipped += 1
                print(f"WARNING [{label}]: line {lineno} has {len(cols)} cols, expected 7 - skipping", file=sys.stderr)
                continue

            qseqid, sseqid, pident, length, evalue_s, bitscore_s, stitle = cols[:7]
            try:
                evalue = float(evalue_s)
                bitscore = float(bitscore_s)
            except ValueError:
                skipped += 1
                print(f"WARNING [{label}]: line {lineno} has non-numeric evalue/bitscore - skipping", file=sys.stderr)
                continue

            candidate = (evalue, -bitscore, sseqid, stitle)
            current = best.get(qseqid)
            if current is None or candidate[:2] < current[:2]:
                best[qseqid] = candidate

    tier_counts = {}
    with open(out_path, "w") as fout:
        fout.write("qseqid\tevalue\tbitscore\ttier\tparsed_name\traw_stitle\n")
        for qseqid, (evalue, neg_bitscore, sseqid, stitle) in sorted(best.items()):
            bitscore = -neg_bitscore
            tier = tier_for_evalue(evalue)
            name = parse_name(stitle)
            tier_counts[tier] = tier_counts.get(tier, 0) + 1
            fout.write(f"{qseqid}\t{evalue}\t{bitscore}\t{tier}\t{name}\t{stitle}\n")

    print(f"[{label}] input lines: {total_lines}, skipped: {skipped}, unique queries: {len(best)}", file=sys.stderr)
    print(f"[{label}] tier breakdown: {tier_counts}", file=sys.stderr)
    print(f"[{label}] wrote {out_path}", file=sys.stderr)


def main():
    process(
        in_path="../results/known_gene_blastx_all_hits.tsv",
        out_path="../results/known_gene_blastx_best_hits_parsed.tsv",
        label="known",
    )


if __name__ == "__main__":
    main()
