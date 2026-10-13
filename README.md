# Known_Genes_Original_14720: per-gene SwissProt re-annotation

Culex tarsalis CtarK1 updated GTF. Re-annotates the 14,720 known genes by running BLASTX on each gene's **own** genomic sequence against SwissProt, choosing the best evidence per gene, and rebuilding the notes in the final GTF.

Last updated: 2026-10-13. Status: pipeline complete; see "Remaining TODOs".
Repo: github.com/yates13/Known_Genes_Original_14720.

---

## Why this pipeline exists

The first annotation run (archived in `results/old_oct6/`) picked one **representative transcript per gene** (`gene_to_representative_transcript.tsv`), ran BLASTX on the transcript, and copied the best SwissProt hit onto the gene.

Problem: many genes shared the same representative, and a long transcript can span several genes. One hit was then copied onto all of them.

Worked example:
- `gene5024` to `gene5028` all had the representative `NODE_6006_length_6914_cov_51.114018_g2800_i0` (6.9 kb).
- BLASTX on that transcript hit Axin, and `note "Axin [blastx_swissprot: extremely_strong]"` was written on all five genes.
- gene5028 (`tig00003883:68689-69185`, minus strand, 497 bp) is actually **RPL32** (60S ribosomal protein L32 / RP49). Its own transcript (`NODE_202888_..._g78408_i2`) hits *D. melanogaster* P04359 at 88% identity over 128 aa, e = 3e-75.

Fix: BLASTX each gene's own sequence (`known_genes.fa`). After the fix gene5024-5028 each carry their own label (NompC-type "Neither inactivation nor afterpotential protein G", SERP2, uS12, RPB7) and gene5028 is "Large ribosomal subunit protein eL32".

Related: GMAP aligned `NODE_176198_..._g78408_i0` across about 46 kb and six genes (gene5023-5028), a spurious long-intron alignment. It is not used by this pipeline.

---

## Environment

- Cluster: CSU Alpine, SLURM (`acpu` partition, `cpu-normal` QoS)
- BLAST+ 2.17.0: `/scratch/alpine/c832500103@colostate.edu/blast_env`
- samtools: `.../conda_envs/hisat2_env/bin/samtools`; bedtools: `.../conda_envs/cellSquito/bin/bedtools` (needs `LD_LIBRARY_PATH=.../hisat2_env/lib`)
- Run every script from `scripts/` (they use `../` paths)

## Inputs

| Input | Source |
|---|---|
| `Spades_Pipeline_Organized/results/basefeatures_realgenes.gtf` (14,720 gene lines; identical copy at `data/basefeatures_realgenes.gtf`) | Earlier pipeline (`Spades_Pipeline_Organized`) |
| `Spades_Pipeline_Organized/data/genome.fa` | CtarK1 genome |
| `results/blastx_db/swissprot.*` | Pre-formatted NCBI SwissProt BLAST db (taxdb files alongside), unpacked from `Spades_Pipeline_Organized/data/swissprot.tar.gz`; no `makeblastdb` was run. **Release: TODO, paste `blastdbcmd -db results/blastx_db/swissprot -info` (date, sequences, letters)** |

## Directory layout

```
Known_Genes_Original_14720/
  data/
    chunks/known_genes.01.fa ... .20.fa    736 genes each, 14,720 total, 0 duplicate IDs
    basefeatures_realgenes.gtf             input GTF (14,720 gene lines)
  logs/                                    SLURM stdout/err (kg_blastx_<jobid>_<task>.out/.err)
  results/
    known_genes.bed, known_genes.fa(.fai)  per-gene coordinates and sequences
    blastx_db/                             SwissProt BLAST db
    blastx_chunks/                         per-gene BLASTX output (one .tsv per chunk)
    old_oct6/                              archived final tables/GTF of the shared-representative run
    old_shared_rep_run/                    old run's inputs and intermediates (representative map,
                                           chunks, BLASTX parts, ID lists); moved 2026-10-13; do NOT use as input
    (Oct 7 no-hit rescue files)            nohit_*, candidate_*, diamond_db/, diamond_test/,
                                           mosquito_pilot/, mosquito_proteins/ (06.x work, rerun on new no-hit list)
  scripts/
    old_scripts_superseded/                replaced scripts and .bak copies
```

---

## Pipeline order

| Step | Script | What it does |
|---|---|---|
| 04.0 | `04.0_extract_known_gene_sequences.sh` | GTF gene lines to BED6, check contigs, `bedtools getfasta -s -nameOnly`, strip `(+)/(-)` from headers. Writes `results/known_genes.bed` and `.fa`. Not rerun (BED verified byte-identical to the original) |
| 04.1 | `04.1_split_genes.sh` | Split `known_genes.fa` into 20 chunks of 736. **Deletes existing chunks; do not rerun while jobs run** |
| 04.2 | `04.2_blastx_array.sh` | SLURM array 1-20, BLASTX each chunk vs SwissProt |
| 04.3 | `04.3_merge_known_gene_blastx.sh` | Merge `results/blastx_chunks/*.blastx.tsv` into `known_gene_blastx_all_hits.tsv` (requires exactly 20 files) |
| 04.4 | `04.4_parse_blastx_hits.py` | Best hit per gene (lowest e-value, tie to highest bitscore), parse name from stitle, assign tier |
| 04.5 | `04.5_join_hits_to_genes.py` | Rewritten for per-gene queries (qseqid = gene_id); output columns unchanged so 05.x runs as before |
| 05.1 | `05.1_build_annotation_evidence.py` | Evidence table, one row per gene (accession, identity, alignment length) |
| 05.1b | `05.1b_choose_best_evidence.py` | One evidence row per gene from the new run and the old transcript run (rules below) |
| 05.2 | `05.2_build_gtf_with_annotation.py` | Write the annotated GTF; reads `annotation_evidence.SwissProt_BLASTX.best.tsv` |

Superseded (old shared-representative approach, in `scripts/old_scripts_superseded/`): `04.1_build_known_gene_representatives.py`, `04.2_prep_known_gene_blastx.sh`, `04.3_known_gene_blastx_array.sh`, plus `.bak` copies of the patched `04.3` merge, old `04.5` join and `05.2`.

### BLASTX settings (04.2)

```
blastx -query known_genes.NN.fa -db ../results/blastx_db/swissprot -evalue 1e-5 \
  -max_target_seqs 5 -max_intron_length 20000 -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle"
```

SLURM: 4 CPUs, 16 GB, 12 h per task. Tasks took 20-33 min each. Output goes to `.tsv.tmp` and is renamed on success, so a failed task leaves no partial `.tsv`.

### Running

```bash
cd scripts
mkdir -p ../logs ../results/blastx_chunks
bash 04.1_split_genes.sh
sbatch --array=1 04.2_blastx_array.sh         # test one chunk first
sbatch --array=2-20 04.2_blastx_array.sh      # then the rest

# after all 20 tasks are COMPLETED:
bash 04.3_merge_known_gene_blastx.sh
python3 04.4_parse_blastx_hits.py
python3 04.5_join_hits_to_genes.py
python3 05.1_build_annotation_evidence.py
python3 05.1b_choose_best_evidence.py
python3 05.2_build_gtf_with_annotation.py
```

Rerun only failed tasks: `sbatch --array=7,12 04.2_blastx_array.sh`.

### Verification after the array

```bash
ls ../results/blastx_chunks/*.blastx.tsv | wc -l          # expect 20
ls ../results/blastx_chunks/*.tmp 2>/dev/null             # expect nothing
sacct -S 2026-10-08 --format=JobID,JobName,State,Elapsed,ExitCode | grep kg_blastx   # all COMPLETED
cat ../results/blastx_chunks/*.blastx.tsv | cut -f1 | sort -u | wc -l   # genes with a hit
```

### Confidence tiers (04.4)

| Tier | E-value | Genes (new run) |
|---|---|---|
| extremely_strong | < 1e-40 | 5,716 |
| very_strong | < 1e-20 | 2,309 |
| strong | < 1e-10 | 1,496 |
| any_hit | < 1e-5 | 918 |

### Evidence selection rules (05.1b)

Gene-span queries include introns, so the new run can give a shorter alignment than a spliced transcript. The old transcript hit is used only where that transcript belonged to the gene alone (`genes_sharing_transcript == no`); borrowed shared-transcript hits are never used.

1. Only one run has a hit: use it (`new_only`, `old_fallback`).
2. Same name: use the new row.
3. New name looks like a transposon and the old one does not: use the old; the reverse: use the new.
4. Otherwise the higher bitscore wins; a tie goes to the new run.

Transposon filter: `transpos|retrotranspos|reverse transcriptase|RNA-directed DNA polymerase|mobile element|jockey|gag-pol|\bpol protein`.

---

## Results (2026-10)

- 14,720 genes queried; 10,439 (70.9%) had a hit in the new run, 4,281 none.
- Compared with the Oct 6 run: 4,104 same name, 1,872 changed, 4,463 gained, 677 lost, 3,604 no hit in either.
- Of the 1,066 changed-name genes whose old transcript was unshared, the old hit had the higher bitscore in 653 and the new hit was equal or higher in 413. Changes among genes with a shared old transcript: 806 (these are the corrections); lost hits: 317 shared, 360 unshared.
- Final GTF: **10,799 genes annotated, 3,921 `annotation_source "none"`**.
- Rule counts: new_only 5,799; same_name 3,574; old_higher_bitscore 611; new_ge_bitscore 338; old_fallback 360; new_chosen_old_is_TE 61; old_chosen_new_is_TE 56.
- Label source in the GTF: 9,772 from the new run, 1,027 from the old transcript run.
- TE-like names in the final set: 193 (217 in the new run alone, 208 in the old run).
- Output GTF has 14,720 lines, same as the input; only gene-line attributes changed.

## Key files

| File | Meaning |
|---|---|
| `results/basefeatures_realgenes.annotated.gtf` | FINAL annotated GTF |
| `results/annotation_lookup.tsv` | gene, source, name, tier |
| `results/annotation_evidence.SwissProt_BLASTX.best.tsv` | chosen evidence, one row per gene |
| `results/known_genes.fa` / `.bed` | per-gene sequences and coordinates used for BLASTX |
| `results/blastx_chunks/known_genes.NN.blastx.tsv` | per-gene BLASTX hits |
| `results/old_oct6/` | archived shared-representative outputs |
| `results/basefeatures_updated_high_confidence.with_notes.gtf` | OLD GTF with the Axin error; do not use |

## Validation done

- Pilot: gene5023-5028 extracted to `pilot_cluster.fa` and run through BLASTX before the full run; Axin did not reappear.
- ID checks: FASTA and GTF both have 14,720 unique gene IDs with none missing; hit IDs not in the FASTA = 0; evidence genes missing from the GTF = 0; the input GTF has no pre-existing `annotation_` or `note` attributes.
- gene5024-5028 labels checked in the final GTF; line count of the annotated GTF equals the input.
- Old vs new comparison (same/changed/gained/lost, shared vs unshared), bitscore comparison and transposon-label check, as reported above.

---

## Not part of the fix

`06.1_diamond_validation.py`, `06.2_candidate_transcripts.py`, `06.4_mosquito_blastx_pilot.sh`: attempts to annotate no-hit genes with DIAMOND and mosquito proteins. Their inputs come from the old no-hit list (Oct 7). 05.1/05.2 do not read any 06.x output.

## Remaining TODOs

- [ ] SwissProt release details for the Inputs table (`blastdbcmd -db results/blastx_db/swissprot -info`).
- [ ] `06.3` is absent from the numbering; check `old_scripts/` and `old_scripts_superseded/`, or renumber `06.4`.
- [ ] Spot-check about 20 of the 611 `old_higher_bitscore` choices.
- [ ] Decide whether to rerun the 06.x rescue on the 3,921 unlabeled genes.
- [ ] Save the validation one-liners as `05.3_validate_annotation.sh`.
- [ ] Update stale step numbers in comments of `04.3`, `04.4`, `04.5`.
