
#!/usr/bin/env python3

"""06.6 - one mosquito-protein (DIAMOND) evidence row per previously unlabeled gene.

Takes every hit <= 1e-10 from Cquinq and Aedes, drops vague names (hypothetical/uncharacterized),

and keeps the highest-bitscore usable hit (ties go to Cquinq). Genes with hits but no usable name

are listed in homology_only_no_usable_name.tsv and get no label. Same columns as the SwissProt

evidence file so 05.2 can read it. Run from scripts/."""

import re, sys

from collections import Counter



DBS = ["Cquinq", "Aedes"]

IN = "../results/mosquito_diamond/{}.tsv"

OUT = "../results/annotation_evidence.Mosquito_DIAMOND.tsv"

HOMOLOGY_ONLY = "../results/mosquito_diamond/homology_only_no_usable_name.tsv"

SOURCE, METHOD = "Mosquito_RefSeq", "DIAMOND"

COLS = ["gene_id", "source", "method", "accession", "name", "evalue", "bitscore",

        "pident", "aln_len", "tier", "transcript_id", "genes_sharing_transcript"]

VAGUE = re.compile(r"hypothetical|uncharacteri[sz]ed|^predicted protein$|^putative protein$|^protein$|^unnamed", re.I)



def clean_name(stitle, sseqid):

    t = stitle.strip()

    if t.startswith(sseqid):

        t = t[len(sseqid):].strip()

    t = re.sub(r"\s*\[[^\[\]]+\]\s*$", "", t)

    t = re.sub(r"^PREDICTED:\s*", "", t, flags=re.I)

    t = re.sub(r"\s+isoform\s+\S+$", "", t, flags=re.I)

    return t.strip()



def accession(sseqid):

    p = sseqid.split("|")

    return p[1] if len(p) >= 3 else sseqid



def tier(ev):

    return "extremely_strong" if ev < 1e-40 else "very_strong" if ev < 1e-20 else "strong" if ev < 1e-10 else "any_hit"



hits = {}   # gene -> list of (bitscore, -rank, row fields)

for rank, db in enumerate(DBS):

    with open(IN.format(db)) as fh:

        for line in fh:

            f = line.rstrip("\n").split("\t")

            if len(f) < 7:

                continue

            g, sseqid, pident, length, ev, bits, stitle = f[:7]

            hits.setdefault(g, []).append((float(bits), -rank, sseqid, pident, length, float(ev), stitle, db))



stats = Counter()

with open(OUT, "w") as out, open(HOMOLOGY_ONLY, "w") as ho:

    out.write("\t".join(COLS) + "\n")

    ho.write("gene_id\tbest_name\tevalue\tbitscore\tdb\n")

    for g in sorted(hits):

        usable = []

        for h in hits[g]:

            name = clean_name(h[6], h[2])

            if name and not VAGUE.search(name):

                usable.append((h, name))

        if usable:

            (bits, _, sseqid, pident, length, ev, _, db), name = max(usable, key=lambda x: (x[0][0], x[0][1]))

            out.write("\t".join([g, SOURCE, METHOD, accession(sseqid), name, f"{ev:g}", f"{bits:g}",

                                 pident, length, tier(ev), g, "no"]) + "\n")

            stats["labeled_" + db] += 1

        else:

            h = max(hits[g], key=lambda x: (x[0], x[1]))

            ho.write(f"{g}\t{clean_name(h[6], h[2])}\t{h[5]:g}\t{h[0]:g}\t{h[7]}\n")

            stats["homology_only"] += 1



print(f"genes with any mosquito hit: {len(hits)}", file=sys.stderr)

for k, v in stats.most_common():

    print(f"  {k}: {v}", file=sys.stderr)

print(f"wrote {OUT}", file=sys.stderr)

