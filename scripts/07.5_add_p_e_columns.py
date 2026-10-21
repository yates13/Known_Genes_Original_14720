#!/usr/bin/env python3
# 07.5_add_p_e_columns.py  (run from scripts/)
import re

BASE = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.gtf"
PV   = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.with_pvalues.gtf"
OUT_P  = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.col10_pvalue.gtf"
OUT_PE = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.col10_pvalue_col11_evalue.gtf"

def attr(s, k):
    m = re.search(rf'(?<![\w]){k} "([^"]*)"', s)
    return m.group(1) if m else "NA"

# gene_id -> label_pvalue, taken from the with_pvalues GTF
pval = {}
with open(PV) as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 9 and f[2] == "gene":
            pval[attr(f[8], "gene_id")] = attr(f[8], "label_pvalue")

n = 0
with open(BASE) as fin, open(OUT_P, "w") as fp, open(OUT_PE, "w") as fpe:
    for line in fin:
        if line.startswith("#"):
            fp.write(line); fpe.write(line); continue
        f = line.rstrip("\n").split("\t")
        gid = attr(f[8], "gene_id") if len(f) >= 9 else ""
        p = pval.get(gid, "NA")
        e = attr(f[8], "annotation_evalue") if len(f) >= 9 else "NA"
        fp.write("\t".join(f + [p]) + "\n")
        fpe.write("\t".join(f + [p, e]) + "\n")
        n += 1
print(n, "lines written")
