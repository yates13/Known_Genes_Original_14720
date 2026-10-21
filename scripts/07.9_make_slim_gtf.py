
#!/usr/bin/env python3

# 07.9_make_slim_gtf.py (run from scripts/)

# Slim 9-column GTF: gene lines keep only the attributes in KEEP (in that order).

import re



IN  = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.gtf"

OUT = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.slim.gtf"

KEEP = ["gene_id", "ID", "Name", "annotation_source", "annotation_method",

        "annotation_accession", "annotation_name", "annotation_evalue",

        "annotation_confidence"]



def attr(s, k):

    m = re.search(r'(?<!\w)%s "([^"]*)"' % k, s)

    return m.group(1) if m else None



n = 0

with open(IN) as fin, open(OUT, "w") as fout:

    for line in fin:

        f = line.rstrip("\n").split("\t")

        if line.startswith("#") or len(f) < 9 or f[2] != "gene":

            fout.write(line); continue

        parts = []

        for k in KEEP:

            v = attr(f[8], k)

            if v is not None:

                parts.append('%s "%s";' % (k, v))

        f[8] = " ".join(parts)

        fout.write("\t".join(f) + "\n")

        n += 1

print(n, "gene lines written to", OUT)

