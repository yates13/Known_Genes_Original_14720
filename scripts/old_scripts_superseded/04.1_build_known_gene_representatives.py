#!/usr/bin/env python3

"""
Step 1 (known genes): pick one representative transcript per known gene_id
from transcripts_confirmed_genes.tsv, so we can blastx those the same way
the novel loci representatives were blasted.

Input:  transcripts_confirmed_genes.tsv
        (bedtools -wao output: transcript BED cols 1-6, gene BED cols 7-12, overlap col 13)
        col4  (idx 3)  = transcript_id  (NODE_..._length_N_...)
        col10 (idx 9)  = gene_id        (e.g. gene1, MT-nbis-gene-2)
        col13 (idx 12) = overlap length in bp

Output:
  representative_transcript_ids_knowngenes.txt
      One transcript ID per line - the chosen representative for each known gene.
  gene_to_representative_transcript.tsv
      Two columns: gene_id <tab> representative_transcript_id
      (keep this - it's the join key back to blastx results later)

Selection rule (mirrors the "longest transcript wins" rule used for novel loci):
  1. Prefer the transcript with the greatest embedded length (from "_length_N_" in its ID)
  2. Tie-break on greatest overlap bp with the gene (column 13)
"""

import re
import sys
from collections import defaultdict

IN_TSV = "transcripts_confirmed_genes.tsv"
OUT_IDS = "representative_transcript_ids_knowngenes.txt"
OUT_MAP = "gene_to_representative_transcript.tsv"

LENGTH_RE = re.compile(r"_length_(\d+)_")
PATH_SUFFIX_RE = re.compile(r"\.path\d+$")

def strip_path_suffix(tid: str) -> str:
    """transcripts_confirmed_genes.tsv IDs carry a .pathN suffix (from the novel
    transcript BED) that transcriptome.fasta headers do NOT have. Strip it so
    downstream fasta extraction / blastx qseqid joins actually match. Same class
    of bug as the .mrna1 vs .path1 mismatch fixed in step 01.9."""
    return PATH_SUFFIX_RE.sub("", tid)

def transcript_length(tid: str) -> int:
    m = LENGTH_RE.search(tid)
    if not m:
        # Shouldn't happen given rnaSPAdes naming, but don't silently misrank if it does
        print(f"WARNING: could not parse length from transcript id: {tid}", file=sys.stderr)
        return -1
    return int(m.group(1))


def main():
    # best[gene_id] = (length, overlap, transcript_id)
    best = {}

    with open(IN_TSV) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 13:
                print(f"WARNING: line {lineno} has {len(cols)} cols, expected 13 - skipping", file=sys.stderr)
                continue

            transcript_id = cols[3]
            gene_id = cols[9]
            try:
                overlap = int(cols[12])
            except ValueError:
                overlap = -1

            if gene_id == "." or gene_id == "-1" or not gene_id:
                # unmatched row (no overlapping gene) - skip
                continue

            length = transcript_length(transcript_id)
            candidate = (length, overlap, transcript_id)

            current = best.get(gene_id)
            if current is None or candidate[:2] > current[:2]:
                best[gene_id] = candidate

    with open(OUT_IDS, "w") as f_ids, open(OUT_MAP, "w") as f_map:
        for gene_id in sorted(best.keys()):
            length, overlap, transcript_id = best[gene_id]
            clean_id = strip_path_suffix(transcript_id)
            f_ids.write(clean_id + "\n")
            f_map.write(f"{gene_id}\t{clean_id}\n")

    print(f"Known genes with a supporting transcript: {len(best)}", file=sys.stderr)
    print(f"Wrote {OUT_IDS} and {OUT_MAP}", file=sys.stderr)

if __name__ == "__main__":
    main()
