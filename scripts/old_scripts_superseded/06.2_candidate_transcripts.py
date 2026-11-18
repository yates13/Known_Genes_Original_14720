
#!/usr/bin/env python3

"""06.2 - candidate transcripts for genes with no SwissProt hit. No searching is done.

Needs only the standard library. Run from scripts/ with the base python."""

import re

import sys

from collections import defaultdict



IN_TSV = "../results/transcripts_confirmed_genes.tsv"

EVIDENCE = "../results/annotation_evidence.SwissProt_BLASTX.tsv"

SEARCHED = "../results/representative_transcript_ids_knowngenes.txt"

HITS = "../results/known_gene_blastx_best_hits_parsed.tsv"

REPS = "../results/gene_to_representative_transcript.tsv"

OUT = "../results/candidate_transcripts_nohit_genes.tsv"

OUT_IDS = "../results/candidate_transcript_ids_unsearched.txt"



LEN_RE = re.compile(r"_length_(\d+)_")

SUFFIX_RE = re.compile(r"\.path\d+$")



annotated = set()

with open(EVIDENCE) as fh:

    next(fh)

    for line in fh:

        annotated.add(line.split("\t")[0])



searched = set()

with open(SEARCHED) as fh:

    for line in fh:

        if line.strip():

            searched.add(line.strip())



hits = {}

with open(HITS) as fh:

    for line in fh:

        f = line.rstrip("\n").split("\t")

        if len(f) < 6:

            continue

        try:

            float(f[1])

        except ValueError:

            continue

        hits[f[0]] = (f[3], f[4])



rep = {}

with open(REPS) as fh:

    for line in fh:

        f = line.rstrip("\n").split("\t")

        if len(f) == 2:

            rep[f[0]] = f[1]



cand = {}

with open(IN_TSV) as fh:

    for line in fh:

        c = line.rstrip("\n").split("\t")

        if len(c) < 13:

            continue

        gene = c[9]

        if gene in (".", "-1", "") or gene in annotated:

            continue

        tx = SUFFIX_RE.sub("", c[3])

        try:

            ov = int(c[12])

            tspan = int(c[2]) - int(c[1])

            glen = int(c[8]) - int(c[7])

        except ValueError:

            continue

        key = (gene, tx)

        if key not in cand or ov > cand[key][0]:

            cand[key] = (ov, tspan, glen)



rows = []

for (gene, tx), (ov, tspan, glen) in cand.items():

    m = LEN_RE.search(tx)

    tlen = int(m.group(1)) if m else -1

    fg = ov / glen if glen > 0 else 0.0

    ft = ov / tspan if tspan > 0 else 0.0

    if tx in hits:

        status = "searched_hit"

    elif tx in searched:

        status = "searched_nohit"

    else:

        status = "unsearched"

    tier, name = hits.get(tx, ("", ""))

    is_rep = "yes" if rep.get(gene) == tx else "no"

    rows.append((gene, tx, tlen, tspan, glen, ov, round(fg, 3), round(ft, 3), is_rep, status, tier, name))

rows.sort()



with open(OUT, "w") as out:

    out.write("gene_id\ttranscript_id\ttx_len\ttx_span\tgene_len\toverlap\tfrac_gene\tfrac_tx\tis_representative\tstatus\tprior_tier\tprior_name\n")

    for r in rows:

        out.write("\t".join(str(x) for x in r) + "\n")



unsearched_ids = sorted({r[1] for r in rows if r[9] == "unsearched"})

with open(OUT_IDS, "w") as out:

    for t in unsearched_ids:

        out.write(t + "\n")





def quant(vals, p):

    vals = sorted(vals)

    return vals[int(p * (len(vals) - 1))] if vals else float("nan")





genes = {r[0] for r in rows}

print(f"no-hit genes with at least one candidate: {len(genes)} (expect 7300)", file=sys.stderr)

print(f"gene-transcript candidate pairs: {len(rows)}", file=sys.stderr)

for s in ("searched_hit", "searched_nohit", "unsearched"):

    sub = [r for r in rows if r[9] == s]

    print(f"{s}: pairs={len(sub)} transcripts={len({r[1] for r in sub})} genes={len({r[0] for r in sub})}", file=sys.stderr)



per_gene = defaultdict(int)

for r in rows:

    per_gene[r[0]] += 1

counts = list(per_gene.values())

print(f"candidates per gene: median={quant(counts, 0.5)} p90={quant(counts, 0.9)} max={max(counts) if counts else 0}", file=sys.stderr)



uns = [r for r in rows if r[9] == "unsearched"]

for label, idx in (("frac_tx", 7), ("frac_gene", 6)):

    qs = [quant([r[idx] for r in uns], p) for p in (0.1, 0.25, 0.5, 0.75, 0.9)]

    print(f"unsearched {label} quantiles (10/25/50/75/90%): {qs}", file=sys.stderr)



print(f"wrote {OUT} and {OUT_IDS} ({len(unsearched_ids)} unsearched transcripts)", file=sys.stderr)

