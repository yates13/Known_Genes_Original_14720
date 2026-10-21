#!/usr/bin/env python3
"""06.8 - species + confidence + coverage attributes on the mosquito-rescue GTF (gene lines only).
Reads the 06.7 rescue GTF and the 06.5b DIAMOND tables (with coverage). Rescue labels are rated by
label bitscore / best bitscore for the gene (>=HIGH high, >=MED medium, else low; DUF/LOC names low).
Low-confidence names become 'uncharacterized protein' (old name and best hit kept as attributes).
'LOW QUALITY PROTEIN:' is stripped from names and flagged. SwissProt names are never changed.
Run from scripts/."""
import re, sys
from collections import Counter

IN_GTF  = "../results/basefeatures_realgenes.annotated.mosquito_rescue.gtf"
OUT_GTF = "../results/basefeatures_realgenes.annotated.mosquito_rescue.v3_confidence.gtf"
COV = "../results/mosquito_diamond_cov/{}.tsv"
DBS = {"Cquinq": "Culex quinquefasciatus", "Aedes": "Aedes aegypti",
       "Agam": "Anopheles gambiae", "Dmel": "Drosophila melanogaster"}
HIGH, MED, SW_HIGH = 0.8, 0.5, 1e-20
VAGUE = re.compile(r"DUF\d|\bLOC\d|uncharacteri[sz]ed|hypothetical", re.I)

def acc(s):
    p = s.split("|")
    return p[1] if len(p) >= 3 else s

def tclean(t, s):
    t = t.strip()
    if t.startswith(s):
        t = t[len(s):].strip()
    t = re.sub(r"\s+OS=.*$", "", t)
    t = re.sub(r"\s*\[[^\[\]]+\]\s*$", "", t)
    return t.replace('"', "'").replace(";", ",").strip()

best, hits = {}, {}   # gene -> (bitscore, acc, title, db); (gene, acc) -> (bitscore, db, scovhsp)
for db in DBS:
    for line in open(COV.format(db)):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        g, s, bs, title, scov = f[0], f[1], float(f[5]), f[6], f[11]
        a = acc(s)
        if g not in best or bs > best[g][0]:
            best[g] = (bs, a, tclean(title, s), db)
        if (g, a) not in hits or bs > hits[(g, a)][0]:
            hits[(g, a)] = (bs, db, scov)

def get(a, k):
    m = re.search(r'(?<!\w)%s "([^"]*)"' % k, a)
    return m.group(1) if m else None

def setv(a, k, v):
    return re.sub(r'(?<!\w)%s "[^"]*"' % k, lambda m: '%s "%s"' % (k, v), a, count=1)

stats = Counter()
with open(IN_GTF) as fin, open(OUT_GTF, "w") as fout:
    for line in fin:
        f = line.rstrip("\n").split("\t")
        if line.startswith("#") or len(f) < 9 or f[2] != "gene":
            fout.write(line); continue
        a = f[8].rstrip()
        if not a.endswith(";"):
            a += ";"
        src, gid, add = get(a, "annotation_source"), get(a, "gene_id"), []
        if src == "SwissProt":
            try:
                conf = "high" if float(get(a, "annotation_evalue")) < SW_HIGH else "medium"
            except (TypeError, ValueError):
                conf = "unverified"
            add.append(("annotation_confidence", conf)); stats[("SwissProt", conf)] += 1
        elif src == "Mosquito_RefSeq":
            accn, name = get(a, "annotation_accession"), get(a, "annotation_name") or ""
            new = re.sub(r"^LOW QUALITY PROTEIN:\s*", "", name)
            h, b = hits.get((gid, accn)), best.get(gid)
            if not h or not b:
                conf, ratio = "unverified", None
            else:
                ratio = h[0] / b[0]
                conf = "low" if (VAGUE.search(new) or ratio < MED) else "medium" if ratio < HIGH else "high"
                add.append(("annotation_species", DBS[h[1]]))
                add.append(("annotation_label_vs_best", f"{ratio:.3f}"))
                add.append(("annotation_protein_coverage", h[2]))
            add.append(("annotation_confidence", conf)); stats[("Mosquito_RefSeq", conf)] += 1
            if new != name:
                add.append(("annotation_flag", "low_quality_model")); stats["lowquality_stripped"] += 1
            if conf == "low":
                add.append(("annotation_name_original", new))
                add.append(("annotation_similar_to", f"{b[1]} {b[2]}"))
                new = "uncharacterized protein"; stats["demoted"] += 1
            if new != name:
                a = setv(a, "annotation_name", new)
                a = re.sub(r'(?<!\w)note "([^"]*)"',
                           lambda m: 'note "%s"' % m.group(1).replace(name, new, 1), a, count=1)
        a += "".join(' %s "%s";' % kv for kv in add)
        f[8] = a
        fout.write("\t".join(f) + "\n")

for k, v in sorted(stats.items(), key=str):
    print(f"{k}: {v}", file=sys.stderr)
print(f"wrote {OUT_GTF}", file=sys.stderr)
