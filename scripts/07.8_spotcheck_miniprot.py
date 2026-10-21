#!/usr/bin/env python3
# 07.8_spotcheck_miniprot.py <ids.txt>   (run from scripts/)
# Spot-check: aligns each gene's label protein (SwissProt or mosquito proteins)
# to the gene sequence with miniprot and reports protein coverage/identity.
# Output: ../results/spotcheck_coverage.tsv (16 columns)
#   status = ok | ok_relaxed | no_alignment | no_label | no_gene_seq | protein_not_found
#
# HOW TO RUN (Alpine)
#   1. Use (base) conda so python3 is 3.14. Do NOT run inside miniprot_env,
#      which has no python and falls back to system python 3.6 (no capture_output).
#        source /scratch/alpine/c832500103@colostate.edu/miniconda3/etc/profile.d/conda.sh
#        conda deactivate        # repeat until the prompt shows (base)
#   2. Put blastdbcmd and miniprot on PATH (not needed if the os.environ line below is kept):
#        export PATH=$PATH:/scratch/alpine/c832500103@colostate.edu/blast_env/bin:/scratch/alpine/c832500103@colostate.edu/conda_envs/miniprot_env/bin
#        which python3 blastdbcmd miniprot
#   3. Check the IDs file before running (never redirect into it with a wrong path;
#      ">" truncates the file even if the command fails):
#        wc -l ../results/spotcheck_ids.txt          # expect 37
#        grep -c '<' ../results/spotcheck_ids.txt    # expect 0 (no placeholders)
#   4. Run:
#        python3 07.8_spotcheck_miniprot.py ../results/spotcheck_ids.txt
#   5. Check the output:
#        wc -l ../results/spotcheck_coverage.tsv                              # expect 38
#        awk -F'\t' '{print NF}' ../results/spotcheck_coverage.tsv | sort | uniq -c   # expect only 16
#
# Miniprot notes: -j2 = insect splice model. Genes with no alignment are retried with
# --outc=0 -k5 -n2 (default --outc=0.1 drops alignments covering <10% of the protein)
# and marked ok_relaxed.
import re, sys, subprocess, tempfile, os

BIN_DIRS = ["/scratch/alpine/c832500103@colostate.edu/blast_env/bin",
            "/scratch/alpine/c832500103@colostate.edu/conda_envs/miniprot_env/bin"]
os.environ["PATH"] += os.pathsep + os.pathsep.join(BIN_DIRS)

GTF = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v2_four_species.col10_pvalue_col11_accession.gtf"
GENES_FA = "../results/known_genes.fa"
SWISS = "../results/blastx_db/swissprot"
MOS = [f"../results/mosquito_proteins/{d}.faa" for d in ["Cquinq", "Aedes", "Agam", "Dmel"]]
OUT = "../results/spotcheck_coverage.tsv"

def read_fa(path, keyfn):
    d, k, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if k: d[k] = "".join(buf)
            k, buf = keyfn(line[1:].strip()), []
        else:
            buf.append(line.strip())
    if k: d[k] = "".join(buf)
    return d

def acc_key(h):
    t = h.split()[0]
    return t.split("|")[1] if "|" in t else t

genes = read_fa(GENES_FA, lambda h: h.split()[0])
mos = {}
for p in MOS:
    for k, v in read_fa(p, acc_key).items():
        mos.setdefault(k, v)

info = {}
for line in open(GTF):
    f = line.rstrip("\n").split("\t")
    if len(f) < 11 or f[2] != "gene": continue
    g = lambda k: (re.search(rf'(?<!\w){k} "([^"]*)"', f[8]) or [None, "NA"])[1]
    info[g("gene_id")] = dict(src=g("annotation_source"), name=g("annotation_name"),
                              acc=f[10], p=f[9])

def swiss_protein(acc):
    r = subprocess.run(["blastdbcmd", "-db", SWISS, "-entry", acc, "-outfmt", "%s"],
                       capture_output=True, text=True)
    return r.stdout.strip().split("\n")[0] if r.returncode == 0 and r.stdout.strip() else None

def protein(acc, src):
    if src == "SwissProt":
        return swiss_protein(acc)
    return mos.get(acc) or mos.get(acc.split(".")[0])

def merged_len(iv):
    iv.sort(); tot, end = 0, -1
    for s, e in iv:
        s = max(s, end)
        if e > s: tot += e - s; end = e
    return tot

ids = [l.strip() for l in open(sys.argv[1]) if l.strip()]
with open(OUT, "w") as out:
    out.write("gene_id\tname\tsource\taccession\tpvalue\tgene_len\tprot_len\tn_aln\tcoverage\tpident_approx\tstatus\n")
    for gid in ids:
        r = info.get(gid)
        if not r or r["acc"] == "NA":
            out.write(f"{gid}\t{r['name'] if r else 'NA'}\t{r['src'] if r else 'NA'}\tNA\t{r['p'] if r else 'NA'}\t\t\t\t\t\tno_label\n"); continue
        gseq, prot = genes.get(gid), protein(r["acc"], r["src"])
        if not gseq or not prot:
            out.write(f"{gid}\t{r['name']}\t{r['src']}\t{r['acc']}\t{r['p']}\t\t\t\t\t\t{'no_gene_seq' if not gseq else 'protein_not_found'}\n"); continue
#        with tempfile.TemporaryDirectory() as td:
#            gf, pf = os.path.join(td, "g.fa"), os.path.join(td, "p.faa")
#            open(gf, "w").write(f">{gid}\n{gseq}\n"); open(pf, "w").write(f">{r['acc']}\n{prot}\n")
#            paf = subprocess.run(["miniprot", "-t2", gf, pf], capture_output=True, text=True).stdout
#        rows = [l.split("\t") for l in paf.splitlines() if l and not l.startswith("##")]
#        if not rows:
#            cov, ident, n = 0.0, "NA", 0
#        else:
#            iv = [(int(x[2]), int(x[3])) for x in rows]
#            cov = merged_len(iv) / len(prot)
#            best = max(rows, key=lambda x: int(x[3]) - int(x[2]))
#            ident, n = f"{int(best[9]) / int(best[10]):.3f}", len(rows)
#        out.write(f"{gid}\t{r['name']}\t{r['src']}\t{r['acc']}\t{r['p']}\t{len(gseq)}\t{len(prot)}\t{n}\t{cov:.3f}\t{ident}\t{'ok' if rows else 'no_alignment'}\n")
        def run_mp(gf, pf, extra=()):
            paf = subprocess.run(["miniprot", "-t2", "-j2", *extra, gf, pf],
                                 capture_output=True, text=True).stdout
            return [l.split("\t") for l in paf.splitlines() if l and not l.startswith("##")]

        with tempfile.TemporaryDirectory() as td:
            gf, pf = os.path.join(td, "g.fa"), os.path.join(td, "p.faa")
            open(gf, "w").write(f">{gid}\n{gseq}\n"); open(pf, "w").write(f">{r['acc']}\n{prot}\n")
            rows, mode = run_mp(gf, pf), "ok"
            if not rows:
                rows, mode = run_mp(gf, pf, ["--outc=0", "-k5", "-n2"]), "ok_relaxed"
        if not rows:
            cov, ident, n, strand, qs, qe, ts, te = 0.0, "NA", 0, "NA", "NA", "NA", "NA", "NA"
            status = "no_alignment"
        else:
            iv = [(int(x[2]), int(x[3])) for x in rows]
            cov = merged_len(iv) / len(prot)
            best = max(rows, key=lambda x: int(x[3]) - int(x[2]))
            ident, n = f"{int(best[9]) / int(best[10]):.3f}", len(rows)
            strand, qs, qe, ts, te = best[4], best[2], best[3], best[7], best[8]
            status = mode
        out.write(f"{gid}\t{r['name']}\t{r['src']}\t{r['acc']}\t{r['p']}\t{len(gseq)}\t{len(prot)}\t{n}\t{cov:.3f}\t{ident}\t{status}\t{strand}\t{qs}\t{qe}\t{ts}\t{te}\n")
print("wrote", OUT)
