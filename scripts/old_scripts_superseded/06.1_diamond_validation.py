#!/usr/bin/env python3
"""06.1 - compare DIAMOND sensitivity modes to the original BLASTX on one chunk (part 001).
No Slurm needed. Run from scripts/ inside the annot_env. Usage: python 06.1_diamond_validation.py [threads]"""
import os
import subprocess
import sys
THREADS = sys.argv[1] if len(sys.argv) > 1 else "4"
Q = "../results/known_gene_chunks/known_gene_representative.part_001.fasta"
BLASTX = "../results/known_gene_blastx_results/blastx_part_001.tsv"
DB = "../results/diamond_db/swissprot"
OUTDIR = "../results/diamond_test"
FMT = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle",
       "qlen", "slen", "qcovhsp", "scovhsp"]
os.makedirs(OUTDIR, exist_ok=True)
def n_queries(path):
    with open(path) as fh:
        return len({line.split("\t")[0] for line in fh if line.strip()})
with open(Q) as fh:
    print("queries in part 001:", sum(1 for l in fh if l.startswith(">")))
print("BLASTX (04.3) queries with a hit:", n_queries(BLASTX))
for mode in ["sensitive", "very-sensitive", "ultra-sensitive"]:
    out = f"{OUTDIR}/part001.{mode}.tsv"
    subprocess.run(["diamond", "blastx", "-q", Q, "-d", DB, "-e", "1e-5", "-k", "5",
                    f"--{mode}", "--threads", THREADS, "--outfmt", "6", *FMT, "-o", out],
                   check=True)
    print(f"DIAMOND {mode} queries with a hit:", n_queries(out), flush=True)
