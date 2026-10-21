# Known_Genes_Original_14720: per-gene SwissProt re-annotation + Mosquito Homology Re-annotation

Culex tarsalis CtarK1 updated GTF. Re-annotates the 14,720 known genes by running BLASTX on each gene's **own** genomic sequence against SwissProt, choosing the best evidence per gene, and rebuilding the notes in the final GTF.

Last updated: 2026-10-21. Status: complete. Final annotation release is the v3_confidence GTF (06.8); 07.x outputs rebuilt on v3; spot-check done. See "Remaining TODOs".
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

Genes lacking SwissProt annotation were searched against four insect
protein databases using DIAMOND:

- Culex quinquefasciatus (Cquinq)
- Aedes spp. (Aedes)
- Anopheles gambiae (Agam)
- Drosophila melanogaster (Dmel)

For each gene, the highest-bitscore non-vague annotation was retained.
Genes annotated only as hypothetical proteins or uncharacterized proteins
were not assigned labels.

Accepted mosquito annotations were merged into the SwissProt-annotated
GTF to create:

results/basefeatures_realgenes.annotated.mosquito_rescue.gtf
---

## Environment

- Cluster: CSU Alpine, SLURM (`acpu` partition, `cpu-normal` QoS)
- BLAST+ 2.17.0: `/scratch/alpine/c832500103@colostate.edu/blast_env`
- samtools: `.../conda_envs/hisat2_env/bin/samtools`; bedtools: `.../conda_envs/cellSquito/bin/bedtools` (needs `LD_LIBRARY_PATH=.../hisat2_env/lib`)
- Python/pandas/numpy: `.../conda_envs/pandas_env` (07.x scripts)
- miniprot: `.../conda_envs/miniprot_env` (spot-check; to be created)
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
    mosquito_diamond/
      Cquinq.tsv
      Aedes.tsv
      Agam.tsv
      Dmel.tsv
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
| 06.1 | `06.1_diamond_validation.py` | Validates candidate no-hit genes and associated transcript evidence before mosquito-homology rescue. |
| 06.2 | `06.2_candidate_transcripts.py` | Identifies candidate transcripts for genes lacking SwissProt annotation; generates transcript lists and rescue inputs. |
| 06.4 | `06.4_mosquito_blastx_pilot.sh` | Pilot DIAMOND/BLASTX search on a subset of unlabeled genes to test mosquito-protein annotation workflows and filtering rules. (No 06.3 exists; see TODOs.) |
| 06.5 | `06.5_mosquito_diamond_unlabeled_genes.sh` | DIAMOND searches of unlabeled genes against Cquinq, Aedes, Agam, Dmel. Per-database hit tables in results/mosquito_diamond/. |
| 06.5a | `06.5a_...` | TODO: one-line description (I don't have this script's contents). |
| 06.5b | `06.5b_...` | Reruns the four DIAMOND searches with query-coverage columns; tables in results/mosquito_diamond_cov/. Best hits identical to 06.5. |
| 06.6 | `06.6_build_mosquito_evidence.py` | Removes vague labels, picks the best name per gene by bitscore, writes annotation_evidence.Mosquito_DIAMOND.tsv. A single pass reproduces all 2,227 labels (including the 1,222 formerly produced by a separate stage-2 filter). |
| 06.7 | `06.7_build_mosquito_rescue_gtf.py` | Adds mosquito labels to genes with annotation_source "none". Output: ...mosquito_rescue.gtf (0 accession differences from the earlier version across all 14,720 genes). |
| 06.8 | `06.8_apply_confidence_rule.py` | Final step. Reads the 06.7 GTF and 06.5b tables; adds annotation_species, annotation_confidence, annotation_label_vs_best, annotation_protein_coverage; demotes low-confidence rescue names; strips "LOW QUALITY PROTEIN:". Writes ...mosquito_rescue.v3_confidence.gtf. Cutoffs (HIGH 0.8, MED 0.5, SW_HIGH 1e-20) are set at the top of the script. Replaces the old command-line annotation_species step. |
| 07.1 | `07.1_transfer_annotation.py` | Copies the per-gene SwissProt annotation (`annotation_evidence.SwissProt_BLASTX.best.tsv`) onto another GTF's gene lines by matching each target gene to a known gene in `known_genes.bed`. Usage: `python3 07.1_transfer_annotation.py <target.gtf> <out.gtf>`. SwissProt evidence only: mosquito-rescue labels are not transferred. Match rules below. |
| 07.2 | `07.2_make_gene_bed_with_evalues.py` | Builds the gene BED with per-database E-values and P-values (see "Gene BED with E-values and P-values"). Run in `pandas_env`. |
| 07.3 | `07.3_pick_spotcheck_genes.py` | Picks 31 genes for protein-alignment spot-checks (every SwissProt and mosquito tier, 5 database disagreements, 5 unlabeled genes with a hit, plus the gene5028 control). Writes `spotcheck_genes.tsv`. |
| 07.4 | `07.4_add_pvalues_to_gtf.py` | Copies a GTF and appends P-value attributes to each gene line (`label_db`, `label_evalue`, `label_pvalue`, `pvalue_<db>`), looked up by gene_id from the 07.2 BED. Non-gene lines unchanged. |
| 07.5 | `07.5_add_p_e_columns.py` | Reads the v2 GTF and the v2 `with_pvalues` GTF; writes two files with extra tab-separated columns after the attributes: P-value in column 10, and a second file with P-value in column 10 and E-value in column 11. See "v2_four_species outputs". |
| 07.6 | `07.6_make_gene_table.py` | Gene tables from the 07.7 GTF. Writes `gene_table.v3_confidence_with_Accession.csv` (18 columns: gene_id through protein_coverage, plus name_original and similar_to) and `gene_table.v3_confidence_short.csv` (14 columns, gene_id through evalue) for easy viewing. Run in `pandas_env`. |
| 07.7 | `07.7_add_accession_column.py` | Adds the accession as an extra column to the col10_pvalue GTF. |
| 07.9 | `07.9_make_slim_gtf.py` | Slim strict 9-column GTF. Same 14,720 gene lines and columns 1-8 as the v3 GTF; column 9 keeps only gene_id, ID, Name, annotation_source, annotation_method, annotation_accession, annotation_name, annotation_evalue and annotation_confidence. Edit `KEEP` at the top of the script to change the attributes. Run from `scripts/`. |

Superseded (old shared-representative approach, in `scripts/old_scripts_superseded/`): `04.1_build_known_gene_representatives.py`, `04.2_prep_known_gene_blastx.sh`, `04.3_known_gene_blastx_array.sh`, plus `.bak` copies of the patched `04.3` merge, old `04.5` join and `05.2`.

**07.1 matching rules** (same strand required; each target gene is compared with `known_genes.bed`):

| Match | Condition | Result |
|---|---|---|
| `id` | same gene_id and overlapping | known gene's SwissProt evidence copied onto the target |
| `overlap` | different id, exactly one known gene lies at least 80% inside the target | evidence copied from that known gene |
| `ambiguous` | target contains 2+ known genes (including same id plus 2+ inside) | no label; `annotation_source "ambiguous_multiple_known_genes"` |
| none | no known gene overlaps | left unchanged, `annotation_source "previous_pipeline"` |

Old `note` and `annotation_*` attributes on matched or ambiguous genes are removed, and the old note is kept as `previous_note`. `gene_id` is never modified, and non-gene lines are copied unchanged. A summary of match counts and the number of known genes claimed by 2+ target genes is printed to stderr.

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

### Confidence rule (06.8)

| Source | Rule | Result |
|---|---|---|
| SwissProt | annotation_evalue < 1e-20 | high (8,283); otherwise medium (2,516). Names never changed. |
| Mosquito_RefSeq | label bitscore / best bitscore for the gene: >= 0.8 high, 0.5-0.8 medium, < 0.5 low | high 1,263 / medium 403 / low 561 |

- Rescue names containing DUF or LOC numbers are low regardless of ratio.
- Low-confidence rescue genes are renamed "uncharacterized protein" (561 genes). The earlier name is kept in `annotation_name_original`, the best hit in `annotation_similar_to`. The `note` attribute is updated to match.
- "LOW QUALITY PROTEIN:" is removed from names (56 genes) and recorded as `annotation_flag "low_quality_model"`.
- Other attributes added: `annotation_species`, `annotation_label_vs_best`, `annotation_protein_coverage`.
- 13,026 genes carry a label (9,546 high, 2,919 medium, 561 low). 12,465 of them have an informative name; the 561 low-confidence genes are "uncharacterized protein". 1,694 genes remain unannotated.
- Cutoffs are provisional and supported by the spot-check below; edit them at the top of the script and rerun to change them.

## Results (2026-10-15)

14,720 genes in the input GTF.

- SwissProt stage: 10,799 genes annotated, 3,921 with no SwissProt annotation.
- Mosquito rescue stage searched previously unlabeled genes against four protein databases: Cquinq, Aedes, Agam, and Dmel.
- DIAMOND rescue evidence: 3,114 genes had at least one mosquito-protein hit; 2,227 received a usable annotation and 887 were classified as homology-only (no usable name).
- Best-hit source among rescued genes: 1,036 Cquinq, 662 Agam, 315 Aedes, and 214 Dmel.
- Updated GTF: 13,026 genes annotated, 1,694 annotation_source "none".
- Annotation coverage improved from 10,799/14,720 (73.4%) after SwissProt annotation alone to 13,026/14,720 (88.5%) after mosquito-homology rescue.
- Mosquito rescue contributed 2,227 additional annotations, reducing the unannotated set by 2,227 genes relative to the SwissProt-only GTF.
- Incorporation of Agam and Dmel recovered 1,005 additional genes beyond the earlier Cquinq+Aedes-only rescue workflow.
- Label sources in the final GTF: 10,799 SwissProt annotations and 2,227 Mosquito_RefSeq annotations.
- Output GTF has 14,720 gene lines, matching the input GTF; only gene-line attributes were modified.
- Final annotation release: results/basefeatures_realgenes.annotated.mosquito_rescue.gtf (a `v2_four_species` version also exists, see below).

## Gene BED with E-values and P-values (07.2)

Output: `results/Cxt_Main_NEW_ANNOT_gene.with_evalues.bed` (14,720 genes plus a header row; 42 tab-separated columns). Built from `basefeatures_realgenes.annotated.mosquito_rescue.gtf`, `annotation_evidence.SwissProt_BLASTX.best.tsv`, `annotation_evidence.Mosquito_DIAMOND.tsv` and `results/mosquito_diamond/*.tsv`. Requested so low-confidence calls can be filtered.

| Columns | Contents |
|---|---|
| 1-6 | BED6: chrom, start, end, gene_id, score, strand |
| 7-8 | annotation_source, annotation_name (from the GTF) |
| 9-16 | `label_*` for the hit that actually named the gene: db, hit, name, evalue (col 12), bitscore, pident, aln_len, tier (col 16). NA for the 1,694 unlabeled genes. A tie across databases is comma-joined in `label_db`. |
| 17-36 | Best hit per database (hit, evalue, bitscore, pident), in the order SwissProt, Cquinq, Aedes, Agam, Dmel (E-values in cols 18, 22, 26, 30, 34) |
| 37-42 | P-values: SwissProt, Cquinq, Aedes, Agam, Dmel, label (`label_pvalue` = col 42) |

- Per-database best hit = lowest E-value, ties to highest bitscore. Agam/Dmel `tr|ACC|...` IDs are stripped to bare accessions to match the evidence tables.
- P-value = 1 - exp(-E), computed as `-expm1(-E)` so tiny values don't underflow. These are per-search values, not corrected for multiple testing across genes.
- Every labeled gene already passes about 1e-5 (largest label E-value is 1e-05). Maximum E-values per database: Cquinq 9.94e-11, Aedes 9.19e-11, Agam 9.52e-06, Dmel 1e-05. Because of this, P is practically equal to E here, and `label_tier`, `label_pident` and `label_aln_len` are the more useful confidence filters.
- NA in the label columns means "not called", not "no similarity": genes with only vague or homology-only hits (about 887) are NA in the label columns but can still have per-database hits in columns 17-36.
- Checks passed: 14,721 lines with 42 columns on every row; label counts match the GTF (10,799 SwissProt, 1,036 Cquinq, 662 Agam, 315 Aedes, 214 Dmel, 1,694 NA); P <= E everywhere; NA in the E and P columns line up; no `nan` strings; gene1 label E = P = 7.62e-49; gene5028 still eL32 (P04359.3).
- Example filter (label E-value of 1e-10 or better): `awk -F'\t' 'NR==1 || ($12!="NA" && $12+0 <= 1e-10)' <bed>`
- Environment: `conda_envs/pandas_env` (pandas, numpy). Run from `scripts/`.

### P-value GTF (07.4)

`results/basefeatures_realgenes.annotated.mosquito_rescue.with_pvalues.gtf`: the rescue GTF with `label_db`, `label_evalue`, `label_pvalue` and `pvalue_<db>` appended to each labeled gene line (genes only carry the databases they have a hit for). Checked: 14,720 lines, non-gene lines identical, 13,026 genes with `label_pvalue`, gene1 and gene5028 correct.

### v3_confidence outputs (07.2, 07.4, 07.5, 07.6, 07.7)

All files are in `results/` and are built from `...mosquito_rescue.v3_confidence.gtf`. The older v2_four_species files are superseded.

| File | Contents |
|---|---|
| `...v3_confidence.gtf` | Final annotated GTF (all 14,720 genes, 9 columns) |
| `...v3_confidence.with_pvalues.gtf` | v3 plus P-value attributes in column 9 (07.4) |
| `...v3_confidence.col10_pvalue.gtf` | Column 10 = label_pvalue (07.5) |
| `...v3_confidence.col10_pvalue_col11_evalue.gtf` | Column 10 = P-value, column 11 = E-value (07.5) |
| `...v3_confidence.col10_pvalue_col11_accession.gtf` | Column 10 = P-value, column 11 = accession (07.7) |
| `Cxt_Main_NEW_ANNOT_gene.v3_confidence.with_evalues.bed` | 42-column gene BED (07.2) |
| `gene_table.v3_confidence_with_Accession.csv` | Gene table, 14,720 rows (07.6) |
| `...v3_confidence.slim.gtf` | Slim 9-column GTF for downstream analysis (07.9); readable by gffread and similar tools |
| `gene_table.v3_confidence_with_Accession.csv` | Full gene table, 14,720 rows, 18 columns (07.6) |
| `gene_table.v3_confidence_short.csv` | Short gene table, 14 columns ending at evalue (07.6) |


Files with 10 or 11 columns are not strict GTF; use the plain v3 or with_pvalues GTF for tools like gffread.

- Unlabeled genes (`annotation_source "none"`) have NA in columns 10 and 11.
- Columns 10 and 11 are identical for every labeled gene in practice, because all calls already pass about 1e-5.
- Checked: 14,720 lines in all three files; 10 (or 11) columns on every row; columns 1-9 identical to v2 (diff = 0); 1,694 NA rows in column 10, equal to the 1,694 `none` genes; no labeled gene has NA in column 11.
- Files with 10 or 11 columns are not strict GTF: tools such as `gffread` expect exactly 9. Use them for filtering and keep the plain v2 or `with_pvalues` GTF for downstream tools.
- Example filter (P-value 1e-10 or better): `awk -F'\t' '$10!="NA" && $10+0 <= 1e-10' <col10 file>`


### Spot-check (07.3, 07.8)

37 genes (31 from `spotcheck_genes.tsv` plus 6 Histone H4 genes) were aligned to their label proteins with miniprot; results in `spotcheck_coverage.tsv`. Coverage (and gene-model length vs protein length) was compared against annotation_confidence.

- SwissProt high: 11 of 11 aligned, including distant vertebrate orthologs (about 0.33-0.40 identity).
- SwissProt medium: 1 of 7 aligned; the rest are short gene models or small/fast-evolving proteins.
- Rescue high: 4 of 6 aligned. Rescue medium: 1 of 2. Rescue low: 1 of 6, and the aligned one (mucin) had 0.26 identity, so the demotions look justified.
- All three Ionotropic receptor genes failed to align at every confidence level, a miniprot limitation on a fast-evolving family, not a labeling error.
- Several low-coverage genes with high identity (gene6563, gene5253) are fragmentary gene models, so protein coverage alone should not be used as a hard filter.

Recommendation: use `annotation_confidence` (high for stringent analyses, high+medium for general use). E-value/P-value cutoffs add little because every label already passes about 1e-5.

## Key files

| File | Meaning |
|---|---|
| `results/basefeatures_realgenes.annotated.mosquito_rescue.gtf` | FINAL annotated GTF (SwissProt + mosquito rescue) |
| `results/basefeatures_realgenes.annotated.mosquito_rescue.v2_four_species.gtf` | v2 version with `annotation_species` and mito/rRNA genes |
| `results/basefeatures_realgenes.annotated.gtf` | SwissProt-only annotated GTF (05.2) |
| `results/...v2_four_species.col10_pvalue.gtf` / `...col10_pvalue_col11_evalue.gtf` | v2 with P-value (col 10) and E-value (col 11) for filtering |
| `results/...with_pvalues.gtf` / `...v2_four_species.with_pvalues.gtf` | GTFs with P-value attributes in column 9 |
| `results/Cxt_Main_NEW_ANNOT_gene.with_evalues.bed` | Gene BED with label, per-database E-values and P-values (42 columns) |
| `results/annotation_lookup.tsv` | gene, source, name, tier |
| `results/annotation_evidence.SwissProt_BLASTX.best.tsv` | chosen evidence, one row per gene |
| `results/annotation_evidence.Mosquito_DIAMOND.tsv` | mosquito rescue evidence |
| `results/known_genes.fa` / `.bed` | per-gene sequences and coordinates used for BLASTX |
| `results/blastx_chunks/known_genes.NN.blastx.tsv` | per-gene BLASTX hits |
| `results/old_oct6/` | archived shared-representative outputs |
| `results/basefeatures_updated_high_confidence.with_notes.gtf` | OLD GTF with the Axin error; do not use |
| `results/...v3_confidence.slim.gtf` | Slim GTF (core attributes + confidence); use this for most downstream work |

## Validation done

- Pilot: gene5023-5028 extracted to `pilot_cluster.fa` and run through BLASTX before the full run; Axin did not reappear.
- ID checks: FASTA and GTF both have 14,720 unique gene IDs with none missing; hit IDs not in the FASTA = 0; evidence genes missing from the GTF = 0; the input GTF has no pre-existing `annotation_` or `note` attributes.
- gene5024-5028 labels checked in the final GTF; line count of the annotated GTF equals the input.
- Old vs new comparison (same/changed/gained/lost, shared vs unshared), bitscore comparison and transposon-label check, as reported above.
- 07.x outputs checked as described in their sections above.

---

## Not part of the fix

`06.1_diamond_validation.py`, `06.2_candidate_transcripts.py`, `06.4_mosquito_blastx_pilot.sh`: early attempts to annotate no-hit genes with DIAMOND and mosquito proteins. Their inputs come from the old no-hit list (Oct 7). 05.1/05.2 do not read any 06.x output; the SwissProt GTF is built without them. The mosquito rescue (06.5-06.7) runs after 05.2 and takes the 05.2 GTF as its input.

## Remaining TODOs

- [ ] SwissProt release details for the Inputs table (`blastdbcmd -db results/blastx_db/swissprot -info`).
- [ ] `06.3` is absent from the numbering; renumber `06.4` or note it.
- [ ] Validation one-liners as `05.3_validate_annotation.sh`.
- [ ] Stale step numbers in comments of `04.3`, `04.4`, `04.5`.
- [ ] `.gitignore`: add `__pycache__/` and `*.bak*`; move `.bak` files into `old_scripts_superseded/`.
