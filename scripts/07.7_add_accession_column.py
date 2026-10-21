#!/usr/bin/env python3
# 07.7_add_accession_column.py  (run from scripts/)
# conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/pandas_env

import re

IN  = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.col10_pvalue.gtf"
OUT = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.col10_pvalue_col11_accession.gtf"

n = 0
with open(IN) as fin, open(OUT, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            fout.write(line); continue
        f = line.rstrip("\n").split("\t")
        m = re.search(r'(?<!\w)annotation_accession "([^"]*)"', f[8])
        f.append(m.group(1) if m else "NA")
        fout.write("\t".join(f) + "\n")
        n += 1
print(n, "lines written")
