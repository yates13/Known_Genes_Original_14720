#!/usr/bin/env python3
# 07.6_make_gene_table.py  (run from scripts/, in pandas_env)
# conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/pandas_env
import re
import pandas as pd

IN  = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.col10_pvalue_col11_accession.gtf"
OUT = "../results/gene_table.v3_confidence_with_Accession.csv"

def attr(s, k):
    m = re.search(rf'(?<!\w){k} "([^"]*)"', s)
    return m.group(1) if m else "NA"

rows = []
with open(IN) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        a = f[8]
        rows.append({
            "gene_id":   attr(a, "gene_id"),
            "contig":    f[0],
            "start":     f[3],
            "end":       f[4],
            "strand":    f[6],
            "source":    attr(a, "annotation_source"),
            "name":      attr(a, "annotation_name"),
            "species":   attr(a, "annotation_species"),
            "accession": attr(a, "annotation_accession"),
            "tier":      attr(a, "annotation_tier"),
            "pident":    attr(a, "annotation_pident"),
            "bitscore":  attr(a, "annotation_bitscore"),
            "Accession":    f[10],
            "pvalue":    f[9],
        })

pd.DataFrame(rows).to_csv(OUT, index=False)
print(len(rows), "rows written to", OUT)
