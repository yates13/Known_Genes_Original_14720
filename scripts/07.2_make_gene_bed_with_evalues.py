#!/usr/bin/env python3
# 07.2_make_gene_bed_with_evalues.py  (run from scripts/)

# Must run this on terminal first!
#conda activate /scratch/alpine/c832500103@colostate.edu/conda_envs/pandas_env
#cd /scratch/alpine/c832500103@colostate.edu/culex_tarsalis_assembly/Known_Genes_Original_14720/scripts
#python3 07.2_make_gene_bed_with_evalues.py

import re, csv
import pandas as pd
import numpy as np

COLS = ["qseqid","sseqid","pident","length","evalue","bitscore","stitle"]
GTF  = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.gtf"
OUT  = "../results/Cxt_Main_NEW_ANNOT_gene.v3_confidence.with_evalues.bed"
DBS  = ["Cquinq","Aedes","Agam","Dmel"]
STRIP = r"^(?:tr|sp)\|([^|]+)\|.*$"   # tr|ACC|ACC_SPECIES -> ACC

def load_db(db):
    d = pd.read_csv(f"../results/mosquito_diamond/{db}.tsv", sep="\t",
                    names=COLS, header=None, quoting=csv.QUOTE_NONE)
    d["db"] = db
    d["sseqid"] = d["sseqid"].str.replace(STRIP, r"\1", regex=True)
    return d

raw = {db: load_db(db) for db in DBS}

# ---- 1. Per-database best hit (lowest evalue, tie -> highest bitscore) ----
def best_per_gene(d, tag):
    d = d.sort_values(["qseqid","evalue","bitscore"], ascending=[True,True,False])
    b = d.drop_duplicates("qseqid").set_index("qseqid")[["sseqid","evalue","bitscore","pident"]]
    b.columns = [f"{tag}_{c}" for c in ["hit","evalue","bitscore","pident"]]
    return b

sw = pd.read_csv("../results/annotation_evidence.SwissProt_BLASTX.best.tsv", sep="\t")
sp = sw.set_index("gene_id")[["accession","evalue","bitscore","pident"]]
sp.columns = [f"SwissProt_{c}" for c in ["hit","evalue","bitscore","pident"]]

per_db = pd.concat([sp] + [best_per_gene(raw[db], db) for db in DBS], axis=1)

# ---- 2. The hit that actually supplied the label ----
sw["label_db"] = "SwissProt"
mq = pd.read_csv("../results/annotation_evidence.Mosquito_DIAMOND.tsv", sep="\t")

# recover which DIAMOND db the accepted hit came from: match (gene, accession, bitscore)
allhits = pd.concat(raw.values())[["qseqid","sseqid","bitscore","db"]]
allhits = allhits.rename(columns={"qseqid":"gene_id","sseqid":"accession"})
allhits["bitscore"] = allhits["bitscore"].astype(float)
mq["bitscore"] = mq["bitscore"].astype(float)
mq = mq.merge(allhits.drop_duplicates(), on=["gene_id","accession","bitscore"], how="left")

# collapse ties (same gene/accession/bitscore found in more than one db)
mq = (mq.groupby(["gene_id","accession"], as_index=False)
        .agg({**{c:"first" for c in mq.columns if c not in ("gene_id","accession","db")},
              "db": lambda s: ",".join(sorted(set(s.dropna())))}))
mq = mq.rename(columns={"db":"label_db"})
mq.loc[mq["label_db"] == "", "label_db"] = "NA"

keep = ["gene_id","label_db","accession","name","evalue","bitscore","pident","aln_len","tier"]
label = pd.concat([sw[keep], mq[keep]]).drop_duplicates("gene_id")
label = label.set_index("gene_id")
label.columns = ["label_db","label_hit","label_name","label_evalue",
                 "label_bitscore","label_pident","label_aln_len","label_tier"]
label["label_aln_len"] = label["label_aln_len"].astype("Int64").astype(object)

# ---- 3. GTF genes -> BED, keep only the label that matches the GTF source ----
rows = []
for line in open(GTF):
    if line.startswith("#"): continue
    f = line.rstrip("\n").split("\t")
    if f[2] != "gene": continue
    g = lambda k: (re.search(rf'{k} "([^"]*)"', f[8]) or [None, "NA"])[1]
    rows.append([f[0], int(f[3])-1, int(f[4]), g("gene_id"), ".", f[6],
                 g("annotation_source"), g("annotation_name")])
bed = pd.DataFrame(rows, columns=["chrom","start","end","gene_id","score","strand",
                                  "annotation_source","annotation_name"])

bed = bed.merge(label, left_on="gene_id", right_index=True, how="left")
# a SwissProt row must not label a gene whose GTF source is Mosquito_RefSeq, and vice versa
is_sp  = bed["annotation_source"].str.contains("SwissProt", na=False)
is_mq  = bed["annotation_source"].str.contains("Mosquito", na=False)
mismatch = (is_sp & (bed["label_db"] != "SwissProt")) | (is_mq & (bed["label_db"] == "SwissProt"))
bed.loc[mismatch | ~(is_sp | is_mq), [c for c in bed.columns if c.startswith("label_")]] = pd.NA
bed = bed.merge(per_db, left_on="gene_id", right_index=True, how="left")
# ---- 4. P-values from E-values: P = 1 - exp(-E), appended as the last columns ----
def evalue_to_p(e):
    return -np.expm1(-pd.to_numeric(e.astype(str), errors="coerce"))   # NA/missing -> NaN

for db in ["SwissProt"] + DBS:
    bed[f"{db}_pvalue"] = evalue_to_p(bed[f"{db}_evalue"])
bed["label_pvalue"] = evalue_to_p(bed["label_evalue"])

bed = bed.fillna("NA")
bed.to_csv(OUT, sep="\t", index=False)
print(len(bed), "genes written")
print(bed["annotation_source"].value_counts())
print(bed["label_db"].value_counts())
