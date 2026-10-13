
#!/usr/bin/env python3

"""05.1b - choose ONE evidence row per gene from the new per-gene BLASTX run and the old

transcript-based run (old rows used only where that transcript belonged to the gene alone).

Rules: same name -> new; TE-like new name loses to non-TE old name; otherwise higher bitscore

wins (ties -> new); a gene with only one source uses it. Provenance is kept in the

transcript_id column (gene id = new run, NODE_... = old run). Run from scripts/."""

import re, sys

from collections import Counter



NEW = "../results/annotation_evidence.SwissProt_BLASTX.tsv"

OLD = "../results/old_oct6/annotation_evidence.SwissProt_BLASTX.tsv"

OUT = "../results/annotation_evidence.SwissProt_BLASTX.best.tsv"

TE = re.compile(r"transpos|retrotranspos|reverse transcriptase|RNA-directed DNA polymerase|"

                r"mobile element|jockey|gag-pol|\bpol protein", re.I)



def load(path, only_unshared=False):

    rows = {}

    with open(path) as fh:

        cols = fh.readline().rstrip("\n").split("\t")

        for line in fh:

            r = dict(zip(cols, line.rstrip("\n").split("\t")))

            if only_unshared and r["genes_sharing_transcript"] != "no":

                continue

            rows.setdefault(r["gene_id"], r)

    return cols, rows



cols, new = load(NEW)

_, old = load(OLD, only_unshared=True)



stats, te_out = Counter(), 0

with open(OUT, "w") as out:

    out.write("\t".join(cols) + "\n")

    for g in sorted(set(new) | set(old)):

        n, o = new.get(g), old.get(g)

        if n and not o:

            pick, why = n, "new_only"

        elif o and not n:

            pick, why = o, "old_fallback"

        elif n["name"] == o["name"]:

            pick, why = n, "same_name"

        else:

            n_te, o_te = bool(TE.search(n["name"])), bool(TE.search(o["name"]))

            if n_te and not o_te:

                pick, why = o, "old_chosen_new_is_TE"

            elif o_te and not n_te:

                pick, why = n, "new_chosen_old_is_TE"

            elif float(o["bitscore"]) > float(n["bitscore"]):

                pick, why = o, "old_higher_bitscore"

            else:

                pick, why = n, "new_ge_bitscore"

        stats[why] += 1

        te_out += bool(TE.search(pick["name"]))

        out.write("\t".join(pick[c] for c in cols) + "\n")



print(f"genes with evidence: {sum(stats.values())}", file=sys.stderr)

for k, v in stats.most_common():

    print(f"  {k}: {v}", file=sys.stderr)

print(f"TE-like names in final set: {te_out}", file=sys.stderr)

print(f"wrote {OUT}", file=sys.stderr)

