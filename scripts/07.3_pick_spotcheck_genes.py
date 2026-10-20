#!/usr/bin/env python3
# 07.3_pick_spotcheck_genes.py  (run from scripts/, pandas_env)
# Must run this on terminal first!
#conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/pandas_env
#cd /scratch/alpine/c832500103@colostate.edu/culex_tarsalis_assembly/Known_Genes_Original_14720/scripts
#python3 07.2_make_gene_bed_with_evalues.py

import pandas as pd
BED = "../results/Cxt_Main_NEW_ANNOT_gene.with_evalues.bed"
OUT = "../results/spotcheck_genes.tsv"
SEED = 13
b = pd.read_csv(BED, sep="\t", na_values="NA")
for c in [c for c in b.columns if c.endswith("_evalue")]:
    b[c] = pd.to_numeric(b[c], errors="coerce")

mos = ["Cquinq","Aedes","Agam","Dmel"]
ev = b[[f"{d}_evalue" for d in mos]]
b["best_mos_db"] = ev.fillna(float("inf")).idxmin(axis=1).str.replace("_evalue","")
b.loc[ev.isna().all(axis=1), "best_mos_db"] = pd.NA
b["group"] = b["label_db"].map(lambda x: "SwissProt" if x == "SwissProt" else "mosquito")

picks = []
# 1. stratified: 3 per tier for SwissProt, 2 per tier for mosquito
lab = b[b["label_db"].notna()]
for (grp, tier), df in lab.groupby(["group","label_tier"]):
    n = 3 if grp == "SwissProt" else 2
    s = df.sample(min(n, len(df)), random_state=SEED).assign(stratum=f"{grp}:{tier}")
    picks.append(s)
# 2. mosquito-labeled genes where another database had a lower evalue than the label's
dis = lab[(lab["group"]=="mosquito") & (lab["best_mos_db"] != lab["label_db"])]
picks.append(dis.sample(min(5, len(dis)), random_state=SEED).assign(stratum="db_disagreement"))
# 3. unlabeled genes that still have a decent per-database hit
un = b[b["label_db"].isna()]
un = un[un[[f"{d}_evalue" for d in mos]].min(axis=1) < 1e-10]
picks.append(un.sample(min(5, len(un)), random_state=SEED).assign(stratum="unlabeled_with_hit"))
# 4. positive control
picks.append(b[b["gene_id"]=="gene5028"].assign(stratum="control"))

out = pd.concat(picks).drop_duplicates("gene_id")
cols = ["gene_id","chrom","start","end","strand","stratum","label_db","label_name",
        "label_hit","label_evalue","label_pident","label_aln_len","label_tier","best_mos_db"]
out[cols].to_csv(OUT, sep="\t", index=False)
print(len(out), "genes ->", OUT)
print(out["stratum"].value_counts())
