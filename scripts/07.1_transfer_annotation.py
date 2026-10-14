
#!/usr/bin/env python3

"""07.1 - copy the per-gene SwissProt annotation (annotation_evidence...best.tsv) onto another GTF.

usage: python3 07.1_transfer_annotation.py <target.gtf> <out.gtf>   (run from scripts/)

Each target gene line is matched to a known gene (known_genes.bed):

  id        same gene_id and overlapping on the same strand

  overlap   different id, exactly one known gene lies >=80% inside the target on the same strand

  ambiguous target contains 2+ known genes (or same id AND 2+ known genes inside): no label

  (none)    no known gene: left unchanged, annotation_source "previous_pipeline"

Old note/annotation_* attributes on matched or ambiguous genes are removed (old note kept as

previous_note). gene_id is never modified. Non-gene lines are copied unchanged."""

import re, sys

from collections import Counter, defaultdict



if len(sys.argv) != 3:

    sys.exit(__doc__)

TARGET, OUT = sys.argv[1], sys.argv[2]

EVID = "../results/annotation_evidence.SwissProt_BLASTX.best.tsv"

KNOWN = "../results/known_genes.bed"

MIN_INSIDE = 0.8



GID = re.compile(r'gene_id "([^"]+)"')

NOTE = re.compile(r'\s*(?<!\w)note "([^"]*)";?')

PREV = re.compile(r'\s*previous_note "[^"]*";?')

ANNOT = re.compile(r'\s*annotation_\w+ "[^"]*";?')



def clean(v):

    return v.replace('"', "'").replace(";", ",").replace("\t", " ").replace("\n", " ").strip()



kidx = defaultdict(list)                      # (contig, strand) -> [(start, end, id)]

with open(KNOWN) as fh:

    for line in fh:

        f = line.rstrip("\n").split("\t")

        kidx[(f[0], f[5])].append((int(f[1]) + 1, int(f[2]), f[3]))



ev = {}

with open(EVID) as fh:

    cols = fh.readline().rstrip("\n").split("\t")

    for line in fh:

        r = dict(zip(cols, line.rstrip("\n").split("\t")))

        ev[r["gene_id"]] = r



def match(gid, c, s, e, st):

    id_hit, inside = None, []

    for ks, ke, kid in kidx.get((c, st), []):

        ov = min(e, ke) - max(s, ks) + 1

        if ov <= 0:

            continue

        if kid == gid:

            id_hit = kid

        if ov / (ke - ks + 1) >= MIN_INSIDE:

            inside.append(kid)

    if id_hit:

        return ("ambiguous", inside) if len(inside) > 1 else ("id", [id_hit])

    if len(inside) == 1:

        return "overlap", inside

    if len(inside) > 1:

        return "ambiguous", inside

    return None, []



def strip_annotation(attrs):

    old = NOTE.search(PREV.sub("", attrs))

    attrs = ANNOT.sub("", NOTE.sub("", PREV.sub("", attrs))).rstrip()

    if not attrs.endswith(";"):

        attrs += ";"

    if old:

        attrs += f' previous_note "{clean(old.group(1))}";'

    return attrs, bool(old)



stats, per_known = Counter(), Counter()

with open(TARGET) as fin, open(OUT, "w") as fout:

    for line in fin:

        c = line.rstrip("\n").split("\t")

        if len(c) < 9 or c[2] != "gene":

            fout.write(line); continue

        m = GID.search(c[8]); gid = m.group(1) if m else ""

        how, kids = match(gid, c[0], int(c[3]), int(c[4]), c[6])

        attrs = c[8].rstrip()

        if how in ("id", "overlap"):

            attrs, had_note = strip_annotation(attrs)

            kid = kids[0]; per_known[kid] += 1

            w = ev.get(kid)

            attrs += f' annotation_match "{how}"; annotation_known_gene_id "{kid}";'

            if w:

                pairs = [("annotation_source", w["source"]), ("annotation_method", w["method"]),

                         ("annotation_accession", w["accession"]), ("annotation_name", w["name"]),

                         ("annotation_evalue", w["evalue"]), ("annotation_bitscore", w["bitscore"]),

                         ("annotation_pident", w["pident"]), ("annotation_tier", w["tier"]),

                         ("annotation_transcript", w["transcript_id"]),

                         ("annotation_shared_transcript", w["genes_sharing_transcript"]),

                         ("note", f'{w["name"]} [{w["method"].lower()}_{w["source"].lower()}: {w["tier"]}]')]

                attrs += "".join(f' {k} "{clean(v)}";' for k, v in pairs)

                stats[f"matched by {how}: annotated"] += 1

            else:

                attrs += ' annotation_source "none";'

                stats[f"matched by {how}: no hit (none)"] += 1

        elif how == "ambiguous":

            attrs, had_note = strip_annotation(attrs)

            attrs += f' annotation_source "ambiguous_multiple_known_genes"; annotation_known_gene_id "{",".join(kids)}";'

            stats["ambiguous (spans 2+ known genes): label withheld"] += 1

        else:

            if not attrs.endswith(";"):

                attrs += ";"

            if "annotation_source" not in attrs:

                attrs += ' annotation_source "previous_pipeline";'

            stats["no known gene (left as is)"] += 1

        c[8] = attrs

        fout.write("\t".join(c) + "\n")



for k, v in sorted(stats.items()):

    print(f"{k}: {v}", file=sys.stderr)

dup = sum(1 for v in per_known.values() if v > 1)

print(f"known genes claimed by 2+ target genes: {dup}", file=sys.stderr)

print(f"wrote {OUT}", file=sys.stderr)

