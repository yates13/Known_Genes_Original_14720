#!/usr/bin/env python3
# 07.4_add_pvalues_to_gtf.py  (run from scripts/, in pandas_env)
# conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/pandas_env
import re
import pandas as pd

GTF = "../results/basefeatures_realgenes.annotated.mosquito_rescue.gtf"
BED = "../results/Cxt_Main_NEW_ANNOT_gene.with_evalues.bed"
OUT = "../results/basefeatures_realgenes.annotated.mosquito_rescue.with_pvalues.gtf"

# BED column -> GTF attribute name
FIELDS = {
    "label_db": "label_db",
    "label_evalue": "label_evalue",
    "label_pvalue": "label_pvalue",
    "SwissProt_pvalue": "pvalue_SwissProt",
    "Cquinq_pvalue": "pvalue_Cquinq",
    "Aedes_pvalue": "pvalue_Aedes",
    "Agam_pvalue": "pvalue_Agam",
    "Dmel_pvalue": "pvalue_Dmel",
}

bed = pd.read_csv(BED, sep="\t", dtype=str, keep_default_na=False).set_index("gene_id")

n = 0
with open(GTF) as fin, open(OUT, "w") as fout:
    for line in fin:
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "gene":
            fout.write(line)
            continue
        gid = re.search(r'gene_id "([^"]*)"', f[8]).group(1)
        attrs = f[8].rstrip()
        if not attrs.endswith(";"):
            attrs += ";"
        if gid in bed.index:
            r = bed.loc[gid]
            for col, key in FIELDS.items():
                if r[col] not in ("NA", ""):
                    attrs += f' {key} "{r[col]}";'
            n += 1
        f[8] = attrs
        fout.write("\t".join(f) + "\n")
print(n, "gene lines updated")
