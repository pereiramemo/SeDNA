# SeDNA — Technical Documentation

**Version documented:** the README cites the release as **v1.0.0 (beta)** via Zenodo
(doi:10.5281/zenodo.23224099). The legacy `old/bin/VERSION.txt` still reports **0.1.0**.
**Code state:** the Nextflow port (`sedna.nf`, `modules/`, `bin/sedna-module-<N>.py`) on `main`,
built on top of commit `f86f2f7` (2026-10-08).
**License:** CC BY-NC-SA 4.0 (`LICENSE`). Developed at NewAtlantis Labs, since acquired by
Ocean BioMetrics.
**Scope:** methods and implementation of the SeDNA bioprospecting pipeline — BGC annotation,
dereplication, protein-domain annotation, BGC embedding representation, and BGC clustering
into Gene Cluster Families (GCFs).

---

## 1. Overview

SeDNA is a pipeline for the bioprospecting analysis of **Biosynthetic Gene Clusters (BGCs)** in
metagenomic data. It takes assembled metagenomic contigs or Metagenome-Assembled Genomes (MAGs)
in FASTA format and produces:

1. a **non-redundant catalog of BGC sequences** (GenBank files + metadata table),
2. a **vector (embedding) representation** of every BGC, and
3. a **clustering of BGCs into GCFs**, i.e. groups of clusters expected to produce similar
   compounds.

The design goal is that BGCs encoding chemically similar products end up close together in the
embedding space, so that (a) grouping them into GCFs is a simple geometric operation, and
(b) distance to reference BGCs (e.g. MIBiG) gives a direct, continuous measure of *novelty*
rather than a binary "hit / no hit".

### 1.1 Pipeline at a glance

```
FASTA (MAGs / assembled contigs), one row per sample in the samplesheet
        |
        |  [1] BGC annotation          MODULE1 antiSMASH | MODULE2 deepBGC | MODULE3 GECCO
        |      three independent, concurrent branches on the same input
        v
   per-tool GBKs + annot_metadata.tsv  (intermediate/module-{1,2,3}/<sample>/)
        |
        |  [2] Dereplication           MODULE4
        |      coordinate-based merging of the three per-tool predictions
        v
   non-redundant BGC catalog          (intermediate/module-4/<sample>/dereplicated/)
        |
        |  [3] Domain annotation       MODULE5  (or MODULE5_PRECOMPUTED, --input_mode from_annot)
        |      pyHMMER hmmsearch vs Pfam-A  ->  cath-resolve-hits
        v
   *_annotdoms_resolved.tsv  (one domain architecture per BGC; intermediate/module-5/<sample>/)
        |
        |  [4] Embedding + [5] clustering   MODULE6 (--embedding_method sif | mean), all samples at once
        |      word2vec domain vectors -> SIF-weighted BGC vector -> BIRCH
        v
   bgc_embeddings-<method>.tsv  +  bgc_clust-<method>-thres<T>.tsv  (bgc_id -> GCF id)
```

Section numbers in this document follow the *methods* ([1]–[5] above). Module numbers
(MODULE1–6) follow the Nextflow processes.

### 1.2 Architecture and execution model

SeDNA is a Nextflow pipeline with the same layout as the sibling MG-Clust pipeline
(`~/workspace/repos/tools/MG-Clust`):

```
sedna.nf                          entry point: --help, --input_mode, --stop_at_module, skip flags
nextflow.config                   params, per-process cpus/memory/maxForks, standard + slurm profiles
modules/sedna-module-<N>.nf       one Nextflow process (or variant) per module
bin/
  sedna-module-<N>.py             all logic of module N; called by the matching .nf
  utils.py                        shared helpers + the former src/bgc_annot and src/bgc_derep code
  download_dbs.py                 pre-fetch Pfam-A / antiSMASH / deepBGC databases into ~/.sedna/db
  download_imgs.py                pre-fetch the six images into the Apptainer cache
  get_bgc_coverage.py             standalone (§9.1)
  query_ai.py                     standalone proof of concept (§9.2)
docker/
  Dockerfile.module-<N>           one micromamba image per module
  dockerbuild_commands.sh
  resources/requirements-module-<N>.yml
resources/
  embeddings/*.pkl                domain embeddings + domain frequency table (staged into MODULE6)
  bgc_class_synonyms.tsv          cross-tool BGC class vocabulary (staged into MODULE4)
tests/data/
  from_fasta.tsv                  samplesheet for the three test MAGs
  oceandna/                       three subsampled OceanDNA MAGs
old/                              the pre-Nextflow version: sedna.sh, install/, original scripts
docs/
```

| Module | Process | Script | Former module(s) | Granularity |
|---|---|---|---|---|
| 1 | `MODULE1` | `bin/sedna-module-1.py` | `run_antismash` | per sample |
| 2 | `MODULE2` | `bin/sedna-module-2.py` | `run_deepbgc` | per sample |
| 3 | `MODULE3` | `bin/sedna-module-3.py` | `run_gecco` | per sample |
| 4 | `MODULE4` | `bin/sedna-module-4.py` | `dereplicate` | per sample |
| 5 | `MODULE5`, `MODULE5_PRECOMPUTED` | `bin/sedna-module-5.py` | `annot_cds`, `import_annot_cds` | per sample |
| 6 | `MODULE6` | `bin/sedna-module-6.py` | `cluster_sif`, `cluster_mean` | once, `.collect()`ed |

Each module runs in its **own container image** (`ghcr.io/pereiramemo/sedna/module-<N>`), built
from its own conda env. This per-module isolation is what makes the "plug-and-play" integration of
third-party BGC predictors practical: antiSMASH 6.1.1 (Python 3.9), deepBGC 0.1.31 (Python 3.7 +
TensorFlow 1.15) and GECCO 0.9.10 (Python 3.12) have mutually incompatible dependency trees and
simply cannot coexist in a single environment.

Modules 1–3 are independent branches of the DAG and run concurrently. The legacy `sedna.sh --module run_all` ran them
in sequence. Their outputs are regrouped per sample with `groupTuple(size: <number of active tools>)`
before module 4, sorted back into tool order and staged as `tool_1`, `tool_2`, …. Module 4 therefore
always merges them in the order antiSMASH → deepBGC → GECCO, where the old `dereplicate.py` used
filesystem walk order.

Usage:

```bash
nextflow run sedna.nf --from_fasta_tsv tests/data/from_fasta.tsv --output_dir results
nextflow run sedna.nf --help
```

`old/bin/sedna.sh` and `old/install/` keep the previous conda-based version runnable for
comparison.

---

## 2. BGC annotation

### 2.1 Which tools, and why

SeDNA does not implement its own BGC detector. It runs **three third-party predictors on the same
input** and merges their output:

| Tool | Version | Detection principle |
|---|---|---|
| [antiSMASH](https://github.com/antismash/antismash) | 6.1.1 | **Rule-based**: curated HMM detection rules for known cluster types |
| [deepBGC](https://github.com/Merck/deepbgc) | 0.1.31 | **Deep learning**: BiLSTM over Pfam-domain sequences, plus a product-class classifier |
| [GECCO](https://github.com/zellerlab/GECCO) | 0.9.10 | **Conditional Random Field** over protein-domain composition |

**Rationale for the ensemble.** The three tools have complementary error profiles:

- antiSMASH is *precise but conservative* — its rules only fire for cluster classes that have
  already been characterised, so genuinely novel chemistry can be missed.
- deepBGC and GECCO are *recall-oriented and rule-free* — they generalise to clusters that match
  no curated rule, which is exactly the regime of interest for marine/environmental metagenomes,
  but they produce more false positives and fuzzier boundaries.

Running them together raises sensitivity (a BGC missed by one tool can be recovered by another)
while the agreement between tools is itself usable as a confidence signal. Because the tool set is
plug-and-play, new predictors can be added later by writing one `run_<tool>.py` wrapper plus one
metadata parser, without touching the downstream modules.

### 2.2 How each tool is invoked

All three wrappers take the same general arguments (`--input_sample`, `--sample_name`,
`--output_dir`, `--nslots`, `--overwrite`) and validate the input FASTA before running.

**antiSMASH** (MODULE1, `bin/sedna-module-1.py`)

```
antismash --cpus <nslots> --databases <db_dir>/antismash --genefinding-tool prodigal-m --taxon bacteria \
          --allow-long-headers --minlength 1000 --minimal \
          --output-dir <out> <input.fasta>
```

- `prodigal-m` (Prodigal metagenomic/anonymous mode) is the default gene finder — appropriate for
  fragmented metagenomic contigs of unknown provenance.
- `--minimal` skips the optional comparative/analysis modules; only detection is needed here.
- Outputs are the `*.region*.gbk` files, moved into `module-1/<sample>/gbks/`.

**deepBGC** (MODULE2, `bin/sedna-module-2.py`)

```
DEEPBGC_DOWNLOADS_DIR=<db_dir>/deepbgc \
deepbgc pipeline --detector deepbgc --classifier product_class \
                 --prodigal-meta-mode --min-nucl 1000 --output <out> <input.fasta>
```

- deepBGC emits a single multi-record GBK, which is split into one file per BGC
  (`bin/utils.py::gbk_splitter`) and renamed to the antiSMASH-style `<contig>.region<NNN>.gbk`
  convention so that all tools produce a uniform file layout.
- Predictions are then filtered **after** the run (`bin/utils.py::gbk_filter`) using
  `--deepbgc_score_thres` (default **0.80**) and `--deepbgc_cds_count_thres` (default **3** CDS).
  Filtering is deliberately done post-hoc rather than via tool arguments, so that relaxing the
  stringency later does not require re-running the tool; discarded clusters are retained under
  `gbks_removed/`.
- The reported BGC class is the `product_class` with the highest score.

**GECCO** (MODULE3, `bin/sedna-module-3.py`)

```
gecco run --genome <input.fasta> --cds 3 --threshold 0.80 --jobs <nslots> --output-dir <out>
```

- `--cds` and `--threshold` come from `--gecco_cds_count_thres` and `--gecco_score_thres`. Unlike
  deepBGC, GECCO filters during the run, not post hoc.

- GECCO's `*_cluster_<n>.gbk` outputs are renamed to `.region<NNN>.gbk`, again for uniformity.

### 2.3 Per-tool metadata

Each wrapper ends by parsing its own GenBank output into a common table,
`module-<N>/<sample>/annot_metadata.tsv`, via a tool-specific parser in
`bin/utils.py` §3 (`antismash_annot_parser`, `deepbgc_annot_parser`,
`gecco_annot_parser`). The schema is identical across tools:

| Field | Meaning |
|---|---|
| `acc` | accession of the source record (contig-derived) |
| `bgc_class` | predicted BGC class/product |
| `start`, `end`, `length` | BGC coordinates on the contig |
| `on_edge` | whether the BGC runs into the contig boundary (fragmented prediction) |
| `contig_id`, `contig_length` | source contig and its length |
| `file` | path to the BGC GenBank file |
| `sample_name`, `tool` | provenance |

This common schema is the contract that makes the dereplication step tool-agnostic.

Notes on normalisation performed by the parsers:

- `on_edge` is provided natively by antiSMASH (`cand_cluster/contig_edge`); for deepBGC it is
  *derived* — a BGC is flagged as on-edge if its start is < 3 bp from the contig start or its end
  is within 30 bp of the contig end. This is why the original assembly FASTA must be passed to the
  parser (it supplies contig lengths).
- `bgc_class` is left as each tool reports it at this stage. Harmonisation across the three
  vocabularies happens later, in the `dereplicate` module — see §3.3.
- Because several BGCs can be predicted on one contig and antiSMASH reuses the contig name as the
  accession, parsers disambiguate rows with a per-contig counter (`acc-<i>`) used as the DataFrame
  index.

---

## 3. BGC dereplication

### 3.1 Why dereplication is needed

The same genomic locus is typically predicted by more than one of the three tools, and the
predicted boundaries rarely agree exactly — one tool may call a 25 kb region where another calls
two adjacent 10 kb regions, or the same cluster with a 2 kb offset. Without dereplication:

- the catalog would contain **2–3 copies of most BGCs**, inflating any count-based diversity or
  abundance estimate;
- the GCF clustering would be dominated by trivially identical duplicates;
- novelty statistics would be biased toward whatever the most prolific tool predicted.

The step therefore produces a **catalog of unique BGC loci** that still integrates evidence from
all tools.

### 3.2 How it works

Dereplication is **coordinate-based**, not sequence-alignment based. It uses only `contig_id`,
`start` and `end` from the per-tool metadata tables — which is what makes it fast and completely
independent of how many tools are integrated.

`bin/sedna-module-4.py` loads the `annot_metadata.tsv` of every directory in `--input_dirs` (in
the order antiSMASH, deepBGC, GECCO; tools that found nothing are skipped) and merges them
**recursively and pairwise**: tables 1 and 2 are merged, the merged catalog becomes the new
table 1, which is then merged with table 3, and so on (`recursive_dereplication`). For each pair:

**(a) Shared contigs** (`bin/utils.py::dereplicate_shared_contigs`).
For every pair of BGCs predicted on the same contig, the relative position of the two intervals
is classified into one of **13 topological cases**, and a representative is selected:

| Case | Topology (tool1 vs tool2) | Decision |
|---|---|---|
| 1 | same start, tool1 shorter | keep tool2 (longer) — *overlapped* |
| 2 | tool1 nested inside tool2 | keep tool2 — *overlapped* |
| 3 | same end, tool1 shorter | keep tool2 — *overlapped* |
| 4 | staggered, tool1 starts later | keep the **longer**; *overlapped* if overlap/longest > threshold, else *partially overlapped* |
| 5 | touching at one base (tool2 then tool1) | keep the longer — *partially overlapped* |
| 6 | disjoint (tool2 before tool1) | keep **both** — *non-overlapped* |
| 7 | same start, tool1 longer | keep tool1 — *overlapped* |
| 8 | tool2 nested inside tool1 | keep tool1 — *overlapped* |
| 9 | same end, tool1 longer | keep tool1 — *overlapped* |
| 10 | staggered, tool1 starts earlier | keep the **longer**; *overlapped* or *partially overlapped* by threshold |
| 11 | touching at one base (tool1 then tool2) | keep the longer — *partially overlapped* |
| 12 | disjoint (tool1 before tool2) | keep **both** — *non-overlapped* |
| 13 | identical coordinates | keep tool1 — *overlapped* |

Two governing rules:

- **Longest-wins.** When two predictions describe the same locus, the longer interval is kept.
  The rationale is that a truncated prediction loses biosynthetic genes, whereas a slightly
  over-extended one only adds flanking context — the former is far more damaging to a
  domain-composition-based representation than the latter.
- **`--overlap_thres` (default 0.75).** For the two staggered cases (4 and 10), the pair counts as
  the *same* BGC only if `overlap / length_of_longest > 0.75`; otherwise the two are recorded as
  *partially overlapped* — the longer one still represents the locus, but the relationship is kept
  in the audit output.

Because the merge is iterative, a BGC that "won" in an earlier comparison can be displaced by a
longer competitor later; the code explicitly removes such superseded entries from the dereplicated
set at each case.

**(b) Non-shared contigs** (`bin/utils.py::dereplicate_non_shared_contigs`). BGCs on contigs seen by only one of the
two tables cannot conflict with anything, so they are carried over directly (case 14).

**(c) Single input table** (`bin/utils.py::sort_single_table`). If only one tool was run, every BGC is
passed straight through (case 14) so that downstream modules receive the same directory layout.

### 3.3 BGC class harmonisation

The three predictors use incompatible product vocabularies, so counting or comparing classes
across tools is meaningless on the raw labels:

| Tool | Vocabulary | Examples |
|---|---|---|
| antiSMASH | ~70 fine-grained rule names | `NRPS`, `NRPS-like`, `T1PKS`, `lanthipeptide-class-i`, `siderophore` |
| GECCO | 8 coarse labels, sometimes compound | `NRP`, `NRP;Polyketide`, `Terpene`, `Unknown` |
| deepBGC | 7 product classes | `Alkaloid`, `NRP`, `Polyketide`, `RiPP`, `Saccharide`, `Terpene`, `Other` |

`bin/utils.py::map_synonyms()` resolves them against `resources/bgc_class_synonyms.tsv` (75 rows;
staged into MODULE4 by Nextflow, overridable with
`--synonyms_file`). The table carries one column per tool plus the target vocabulary:

```
bgc_class_antismash   bgc_class_gecco    bgc_class_deepbgc   bgc_class_formatted
NRPS                  NRP                NRP                 nrp
NRPS-like             NRP                NRP                 nrp
NRPS                  NRP;Polyketide     NRP                 nrp
```

The `tool` field of each BGC selects the lookup column (`bgc_class_<tool>`), and the resulting
value is attached to the metadata row as a new `bgc_class_formatted` field. The step runs after
the GBKs are copied and before the metadata tables are written (`sedna-module-4.py` §3.7),
so `dereplicated/annot_metadata.tsv` carries **both** the tool's original `bgc_class` and the
harmonised `bgc_class_formatted`.

The mapping is applied only to the dereplicated dictionary; §11.9 and §11.10 describe two sharp
edges in the current implementation (an unguarded lookup, and leakage of the new column into the
audit tables).

### 3.4 Output

`bin/sedna-module-4.py` writes four sibling directories under `module-4/<sample>/`:

```
module-4/<sample>/
  dereplicated/           <- the catalog: one entry per unique locus
    annot_metadata.tsv
    gbks/                 <- copies of the winning per-tool GBK files, named <bgc_id>.gbk
  overlapped/             <- audit: predictions that matched another prediction
  partially_overlapped/   <- audit: staggered matches below --overlap_thres
  non_overlapped/         <- audit: predictions with no counterpart
```

Each retained BGC receives a stable, self-describing identifier:

```
bgc_id = <tool>__<sample_name>__<contig_id>__<start>__<end>
```

and `dereplicated/annot_metadata.tsv` gains three columns on top of the per-tool schema of §2.3:
`bgc_id`, `link` (the GBK path relative to the table, `gbks/<bgc_id>.gbk`) and `bgc_class_formatted` (§3.3). The three
audit tables get `bgc_id` and `link` but, by the aliasing described in §11.10, an inconsistent
`bgc_class_formatted`. The provenance of every sequence is kept in the `tool` and `file` columns
and in the `bgc_id` prefix. Before the Nextflow port the catalog held symlinks into the per-tool
directories. Those would point into upstream task work dirs, which are neither mounted into the
MODULE4 container nor kept after `nextflow clean`, so the catalog now holds copies.

---

## 4. BGC domain annotation

Clustering operates on **protein-domain composition**, not on nucleotide or protein sequence
similarity. Each BGC is reduced to the ordered list of Pfam domains encoded in its CDSs — the
"sentence of domains" that the embedding step consumes.

### 4.1 `annot_cds` — pyHMMER + cath-resolve-hits

`bin/sedna-module-5.py` (MODULE5) performs three steps for every `*.gbk` in the dereplicated catalog:

**Step 1 — CDS extraction.** `bin/utils.py::extract_cds` writes the translated CDSs of each BGC to a
FASTA file, using `--cds_id` (default `locus_tag`; `protein_id` for MIBiG-style records) as the
sequence identifier.

**Step 2 — HMM search with [pyHMMER](https://github.com/althonos/pyhmmer)** (0.10.x):

```python
pyhmmer.hmmer.hmmsearch(hmms, seqs, domE=evalue_thres, bit_cutoffs="gathering", cpus=nslots)
```

- Database: **Pfam-A** (`<db_dir>/pfam/Pfam-A.hmm`, downloaded by `bin/download_dbs.py`).
- `--cut_ga` (**default on**) uses Pfam's per-family curated **gathering thresholds** rather than a
  single global E-value — the family-specific cutoffs are what Pfam itself recommends and they
  behave much more consistently across domain families of different information content.
  `--no-cut_ga` switches to a plain E-value cutoff (`--evalue_thres`, default `1e-3`).
- **Why pyHMMER and not the `hmmsearch` binary:** pyHMMER binds HMMER3 directly in-process
  (Cython), so there is no subprocess launch, no writing/parsing of intermediate files per BGC and
  no external HMMER installation to manage. With thousands of small per-BGC searches, the
  per-invocation overhead of the CLI dominates; pyHMMER also parallelises internally over `cpus`.
  Output is written in standard `domtblout` format, so downstream tooling is unchanged.

**Step 3 — domain architecture resolution with
[cath-resolve-hits](https://github.com/UCLOrengoGroup/cath-tools)** (`cath-tools` 0.16.5):

```
cath-resolve-hits --input-format hmmer_domtblout <file>.domtblout > <file>_annotdoms_resolved.tsv
```

A raw `domtblout` contains **mutually overlapping hits on the same protein region**: homologous
Pfam families compete for the same stretch of sequence, and clan members frequently all hit.
Taking those hits at face value would (i) count the same region several times and (ii) insert
family-level noise into the domain "sentence", both of which distort the resulting BGC vector.
`cath-resolve-hits` solves this by computing the **highest-scoring, non-overlapping subset** of
hits along each protein (a dynamic-programming optimisation over hit boundaries and scores) — i.e.
the single best-supported domain architecture per CDS.

Output format (space-separated, `#`-commented header):

```
#FIELDS query-id match-id score boundaries resolved cond-evalue indp-evalue
OceanDNA-b9989_00013|ctg13_10|238-1420|1  Arginosuc_synth   139.1  3-163    3-163    6.9e-45  1.5e-43
OceanDNA-b9989_00013|ctg13_10|238-1420|1  Arginosuc_syn_C   110.9  171-386  171-386  3.6e-36  7.5e-35
```

Column 2 (`match-id`, the Pfam family **name**, e.g. `Condensation`, `PP-binding`) is the only
field the embedding step reads. Files are named `<bgc_id>_annotdoms_resolved.tsv`, one per BGC.

### 4.2 `import_annot_cds` — reusing existing annotations

`MODULE5_PRECOMPUTED` (`--input_mode from_annot`; `bin/sedna-module-5.py --gff --domtblout`) produces the same `*_annotdoms_resolved.tsv` files from annotation that
already exists elsewhere: given the module-4 catalog's `annot_metadata.tsv` (BGC coordinates), a `--gff` of gene features
and a `--domtblout` of a genome-wide hmmsearch, it selects the CDSs falling inside each BGC and
runs the resolution step on that subset. This avoids re-running HMM searches when the host
pipeline has already annotated the MAGs.

---

## 5. BGC embedding representation

### 5.1 The analogy

The representation treats BGCs as text:

| NLP | SeDNA |
|---|---|
| word | Pfam domain |
| sentence | BGC (its ordered list of domains) |
| corpus | all BGCs of a reference database |
| word embedding | domain embedding (100-d) |
| sentence embedding | BGC embedding (100-d) |

The premise is distributional: **domains that recur in the same biosynthetic context encode
related chemistry**, exactly as words sharing contexts share meaning. A BGC is then a point in a
continuous space where Euclidean distance approximates biosynthetic (and, by proxy, chemical)
similarity — as opposed to the all-or-nothing matching of a rule- or BLAST-based comparison.

### 5.2 Training the domain embeddings (word2vec)

The shipped model is `resources/embeddings/dom_embeddings_mtx_mibig_pfam.pkl`
(staged into MODULE6 by Nextflow), a plain `dict` mapping a Pfam family name to a
list of 100 floats. **1,262 keys**, 100 dimensions.

Training (R, `word2vec` package; development scripts live outside this repository, in
`dev/new_atlantis/dev/clustering/scripts/R/2.2-create_embeddings_mibig_pfam_evalue_thres.R`):

- **Corpus:** all BGCs of **MIBiG 3.1**, each annotated with the same `annot_cds` procedure
  (pyHMMER vs Pfam-A, then cath-resolve-hits) and flattened to one sentence of domain names per
  BGC, in gene order and, within a gene, in domain order. Domains are **not** deduplicated during
  training — repetition (e.g. tandem PKS modules) is real signal.
- **Model:**

  ```r
  word2vec(x = bgc2dom_clean_list, type = "skip-gram",
           dim = 100, iter = 50, window = 8, min_count = 5)
  ```

  | Parameter | Value | Reason |
  |---|---|---|
  | `type` | skip-gram | better for infrequent tokens than CBOW; most Pfam domains are rare |
  | `dim` | 100 | compact enough for fast distance computation over large catalogs |
  | `window` | 8 | a context of ±8 domains spans several neighbouring genes, capturing module/operon structure rather than only within-gene adjacency |
  | `iter` | 50 | the corpus is small (thousands of sentences), so many epochs are affordable |
  | `min_count` | 5 | drops domains too rare to estimate a vector for; yields the 1,262-token vocabulary |

- Why **MIBiG** as the training corpus: it is the curated reference set of *experimentally
  characterised* BGCs, so the contexts the model learns are real biosynthetic contexts rather than
  prediction artefacts.

### 5.3 From domain vectors to a BGC vector — SIF (`cluster_sif`)

`bin/sedna-module-6.py --embedding_method sif` (formerly `cluster_sif.py`) implements **Smooth Inverse Frequency** weighting (Arora, Liang & Ma, ICLR
2017, *"A Simple but Tough-to-Beat Baseline for Sentence Embeddings"*,
<https://openreview.net/pdf?id=SyK00v5xx>) — a sentence-embedding baseline that is competitive
with far heavier models while remaining fully deterministic and training-free at inference time.

**(a) Frequency weights.** From `resources/embeddings/domain2counts_mibig3.1_vs_pfam.pkl` — a table
of Pfam domain occurrence counts over the same MIBiG 3.1 corpus (3,397 domains, 91,075
occurrences). Domains occurring fewer than **10** times are dropped from the weight table (656
remain), and for the rest:

```
w(d) = alpha / (alpha + freq(d)),     alpha = 0.001 (--alpha)
```

With the shipped table this spans `w = 0.0088` for the most frequent domain (`PP-binding`,
9,451 occurrences) up to `w = 0.894` for domains at the count-10 boundary. The effect is a soft
down-weighting of biosynthetic "stop words" — domains such as `PP-binding`, `ketoacyl-synt`,
`AMP-binding` or `Condensation` that appear in nearly every NRPS/PKS cluster and therefore carry
almost no discriminative information — while letting the rarer, cluster-defining domains dominate
the vector. This is precisely the failure mode of a plain unweighted mean, where a handful of
ubiquitous domains pull every NRPS/PKS BGC toward the same region of space.

**(b) BGC vector.** For each BGC (`map_domain2embeddings`):

1. deduplicate the domain list (each distinct domain contributes once, so long repetitive
   assembly-line clusters do not swamp the vector),
2. sum `w(d) · v(d)` over the domains that exist in the embedding vocabulary,
3. L2-normalise the sum; a BGC with no mappable domain yields an all-`NaN` vector and is set aside.

**(c) Common-component removal.** `remove_pc` fits a 1-component PCA over the whole set of BGC
vectors and projects it out:

```python
output = embeddings - embeddings.dot(pc.T) * pc
```

This is the second half of the SIF recipe. The first principal component of a set of
frequency-weighted sentence vectors is dominated by whatever is common to *all* of them — here,
"being a BGC at all" — and contributes a large shared offset that suppresses the relative
differences the clustering needs. Removing it sharpens the between-BGC contrast.

The result is written to `bgc_embeddings-sif.tsv` (rows = BGCs, 100 columns) and is directly reusable
for downstream novelty and diversity analyses.

### 5.4 `cluster_mean` — the simpler alternative

`--embedding_method mean` (formerly `cluster_mean.py`) is the earlier/baseline variant: it maps domains to vectors and takes the
**unweighted arithmetic mean**, with no SIF weighting, no L2 normalisation and no principal
component removal, then clusters with BIRCH exactly as `cluster_sif` does. It is retained for
comparison; `cluster_sif` is the recommended module.

---

## 6. BGC clustering into GCFs — BIRCH

### 6.1 Configuration

MODULE6 uses scikit-learn's
[`Birch`](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.Birch.html)
(Balanced Iterative Reducing and Clustering using Hierarchies, Zhang et al., SIGMOD 1996):

```python
birch = Birch(n_clusters=None, compute_labels=False, copy=False)   # sif; mean: Birch(n_clusters=None)
birch.threshold = <--birch_threshold>
birch.branching_factor = <number of BGCs>
birch.fit(X)
labels = birch.predict(X)
```

- `n_clusters=None` — **no global clustering step**. The leaf subclusters of the CF-tree *are* the
  GCFs. There is no need to state in advance how many gene cluster families a metagenome contains
  (which is the actual unknown), and no arbitrary cut of a dendrogram.
- `threshold` (`--birch_threshold`, default 1; 0.5 in the README example) is the **only** parameter that
  matters biologically: it is the maximum radius of a leaf subcluster in embedding space, i.e. the
  granularity of a GCF. Lower = tighter, more numerous families.
- `branching_factor` is set to the number of BGCs, i.e. no artificial limit on the number of
  children per node.
- `predict` then assigns each BGC to the nearest subcluster centroid.

### 6.2 Why BIRCH

| Property | Consequence for SeDNA |
|---|---|
| **Linear time, single pass, O(n) memory** | Metagenomic BGC catalogs reach 10^5–10^6 clusters. Pairwise-distance methods (hierarchical clustering, affinity propagation, the all-vs-all network used by BiG-SCAPE) are O(n²) in time *and* memory and become infeasible at that scale. BIRCH builds a Clustering Feature tree in one pass and never materialises a distance matrix. |
| **No `k` required** | The number of GCFs is the result, not an input — unlike k-means/BiG-SLiCE-style approaches where `k` must be guessed. |
| **One interpretable parameter** | The `threshold` is a radius in embedding space, so it maps directly onto "how similar must two BGCs be to be called the same family", and can be swept to produce a family hierarchy at several resolutions. |
| **Incremental / streamable** | The CF-tree can absorb new points without re-clustering from scratch — the right property for a catalog that grows as new samples are sequenced. |
| **Compact summaries (CF triples)** | Each subcluster is stored as (n, linear sum, squared sum) rather than as its members, which is what keeps memory flat. |
| **Outlier tolerance** | Singletons naturally end up in their own leaf, which is the desired behaviour for a genuinely novel BGC — it should form its own family rather than be forced into an existing one. |

Its main limitation — BIRCH assumes roughly spherical, Euclidean-metric clusters and is sensitive
to insertion order — is acceptable here because the SIF embedding is explicitly built to be a
Euclidean space with L2-normalised vectors.

### 6.3 Handling BGCs without domain annotation

BGCs whose resolved-domain file is empty, or whose domains are all outside the embedding
vocabulary, have no vector. Rather than being dropped, they are appended to the output as
**singleton GCFs** with cluster ids continuing past the maximum assigned by BIRCH. This keeps the
catalog and the clustering in one-to-one correspondence.

### 6.4 Output

`bgc_clust-<method>-thres<T>.tsv` (e.g. `bgc_clust-sif-thres0.5.tsv`, in `intermediate/module-6/`
and symlinked into `workables/`) is a two-column table:

```
bgc_id    cluster_id
```

where `cluster_id` is the GCF identifier. Combined with `annot_metadata.tsv` this supports GCF
composition per sample, biosynthetic diversity estimates, and — when reference BGCs (MIBiG) are
included in the same run — functional transfer and novelty assessment: a GCF containing no
reference member is a candidate novel family.

---

## 7. Embedding fine-tuning by combinatorial optimization

> **Status: described from the author's account; no implementation of this step exists in this
> repository, and it is not applied to the shipped model.** The file
> `resources/embeddings/dom_embeddings_mtx_mibig_pfam.pkl` was verified to be bit-identical to the
> raw word2vec output of the training script (`dom_embeddings_mtx_mibig_pfam.tsv`), so no learned
> re-weighting is currently baked into it, and neither embedding method of `bin/sedna-module-6.py`
> applies per-dimension weights. This section records the method for documentation purposes and
> should be reconciled with the original analysis scripts before being cited as implemented.

### 7.1 Motivation

word2vec is trained on a purely *syntactic* objective — predicting neighbouring domains. Nothing in
that objective forces the resulting space to align with **chemical** similarity of the encoded
products, which is the property the GCF clustering actually needs. The 100 dimensions are also not
equally informative: some capture chemically meaningful variation, others capture corpus artefacts.

### 7.2 Approach

A supervised, post-hoc re-weighting of the embedding space, fitted against known
BGC→compound pairs:

1. **Ground truth.** BGC sequences and their associated compounds were taken from the **MIBiG**
   database, which links each characterised BGC to the structure(s) of its product.
2. **Chemical distance.** For each pair of compounds, a **Tanimoto dissimilarity**
   (`1 − Tanimoto coefficient` over molecular fingerprints) was computed — the standard measure of
   structural dissimilarity in cheminformatics.
3. **Genomic distance.** For each corresponding pair of BGCs, the **Euclidean distance** between
   their embeddings.
4. **Objective.** A weight vector applied to the embedding (per-dimension scaling of the 100
   coordinates, equivalently a diagonal metric on the space) was searched so as to **maximise the
   agreement between the inter-BGC Euclidean distances and the inter-compound Tanimoto
   dissimilarities** — i.e. BGC pairs making similar molecules are pulled together, pairs making
   dissimilar molecules pushed apart.
5. **Search.** The weights were found by **combinatorial optimization** over the weight
   configurations, rather than by gradient descent — the objective is a rank/correlation-style
   criterion computed over all pairs and is not conveniently differentiable.

### 7.3 Open points to confirm against the original scripts

- Whether the weights are **per embedding dimension** (100 weights, a diagonal metric) or
  **per domain** (a learned replacement for the SIF weights). The text above assumes
  per-dimension; the author's recollection is uncertain on this point.
- The exact objective: maximised correlation (Pearson/Spearman/Mantel) between the two distance
  matrices, or a margin-style criterion.
- The fingerprint type and Tanimoto variant used for the compound distances, and which MIBiG
  release / subset of BGC–compound pairs was used.
- The optimizer (e.g. simulated annealing, genetic algorithm, exhaustive/greedy search) and its
  stopping criterion.
- Whether a held-out split was used to check that the fitted weights generalise beyond MIBiG.

---

## 8. Benchmarking

The clustering methodology was evaluated (in the development repository, outside this codebase)
against the two established GCF tools:

- **Reference standard:** a manually curated grouping of MIBiG BGCs into GCFs, taken from the
  supplementary material of the BiG-SCAPE/CORASON publication
  (Navarro-Muñoz et al., *Nat. Chem. Biol.* 2020, doi:10.1038/s41589-019-0400-9), reformatted into
  291 BGCs with domain annotation (313 before filtering).
- **Competitors:** BiG-SCAPE (alignment/network based) and BiG-SLiCE (feature-vector + k-means).
- **Metric:** **V-measure** (harmonic mean of homogeneity and completeness) between each tool's
  clustering and the reference, computed over a sweep of each tool's cutoff parameter — for SeDNA,
  the BIRCH `threshold` (swept from 0.01 to 0.7 in the benchmark scripts).

Sweeping the cutoff rather than reporting a single value is important: the three methods have
non-comparable similarity scales, so only the best-achievable agreement per method is a fair
comparison.

---

## 9. Auxiliary modules

### 9.1 `get_bgc_coverage`

`bin/get_bgc_coverage.py` estimates read coverage (i.e. relative abundance) for each BGC in the
catalog, from either a BAM file of reads mapped to the assembly or a per-contig coverage table
(e.g. produced by VEBA), using `bedtools`/`pysam` over the BGC coordinates in `annot_metadata.tsv`.
This turns presence/absence GCF profiles into quantitative ones for cross-sample comparison.

### 9.2 `query_ai` — proof of concept, not wired in

`bin/query_ai.py` takes a JSON file of BGC descriptors (`domains`, `class`, `mibig_hit`,
`taxonomy`, `environment`), fills them into a fixed prompt template asking for a hypothesis about
the cluster's biotechnological potential, and POSTs it to a chat-completions endpoint over
`urllib` (no SDK dependency, stdlib only). The key comes from `--openai-api-key` or
`$OPENAI_API_KEY`; `--test` runs two assertions on the prompt builder and exits without a network
call.

Status and caveats:

- It is **not a pipeline module**: no process calls it, and it has no image of its own (it is
  stdlib only).
- It consumes a hand-written JSON file; nothing in the pipeline emits that file, so there is no
  path yet from `annot_metadata.tsv` + `bgc_clust-*.tsv` to a `query_ai` input.
- The default model is `gpt-4` against `api.openai.com`. Both the provider and that model id are
  hard-coded defaults from the original POC and worth revisiting before this becomes a supported
  module.
- Output is free text printed to stdout; it is not captured, not attached to any BGC record, and
  not reproducible run to run.

---

## 10. Installation and data dependencies

```bash
bash docker/dockerbuild_commands.sh   # builds the six module images (run from the repo root)
bin/download_imgs.py                  # or: cache the images for Apptainer (slurm profile)
bin/download_dbs.py                   # fetches Pfam-A, antiSMASH and deepBGC databases
```

`bin/download_dbs.py` populates `--db_dir` (default `~/.sedna/db`, which must match `params.db_dir`).
The directory is bind-mounted into MODULE1, MODULE2 and MODULE5*; nothing is baked into the images:

| Path | Content | Used by |
|---|---|---|
| `<db_dir>/pfam/Pfam-A.hmm` | Pfam-A HMM library (EBI, current release), fetched on the host | MODULE5 |
| `<db_dir>/antismash/` | antiSMASH databases (`download-antismash-databases`, run in the module-1 image) | MODULE1 (`--databases`) |
| `<db_dir>/deepbgc/data/` | deepBGC models (`deepbgc download`, run in the module-2 image) | MODULE2 (`DEEPBGC_DOWNLOADS_DIR`) |

The layout matches the legacy `~/.local/share`, so `--db_dir ~/.local/share` reuses databases fetched
by `old/install/get_databases.sh`. The embedding pickles and `bgc_class_synonyms.tsv` stay in
`resources/` and are staged into MODULE6 and MODULE4 as ordinary Nextflow inputs.

The previous conda-based installer is kept in `old/install/`.

---

## 11. Implementation notes and caveats

These are behaviours of the current code worth being aware of when interpreting results or
modifying the pipeline.

1. **SIF weight fallback inverts the SIF intent for rare domains.** In
   `bin/sedna-module-6.py::map_domain2embeddings` (formerly `cluster_sif.py`), a domain absent from the filtered weight table receives
   `weights.get(domain, 1e-3)`. The weight table only keeps domains with count ≥ 10 (656 of 3,397),
   so of the 1,262 domains in the embedding vocabulary, **580 rare domains (plus 26 keys with no
   count entry, including word2vec's `</s>` token) receive the 1e-3 fallback** — *lower* than the
   weight of the most frequent domain in the corpus (0.0088). Under SIF, rare domains should
   receive the *highest* weights (up to 0.894 here). The practical effect is that rare — and
   therefore potentially most diagnostic — domains are nearly excluded from the BGC vector.
2. *(Docstring fixed in the Nextflow port.)* **Sum vs. weighted mean.** The old `cluster_sif.py`
   said "weighted mean", but it sums the weighted domain vectors and divides by the L2 norm, and its
   accumulated `total_weight` variable was unused. This is consistent with the SIF paper up to a
   scale factor, and the L2 normalisation makes the difference immaterial. The port removed the
   unused variable and the docstring now says "weighted sum"; the computation is unchanged.
3. **Embedding corpus vs. default annotation settings.** The shipped embedding matrix and the
   frequency table were derived from the MIBiG 3.1 Pfam annotation produced with an **E-value**
   threshold, while MODULE5 (formerly `annot_cds.py`) defaults to **gathering thresholds** (`--cut_ga`). The two
   regimes produce slightly different domain vocabularies (26 embedded domains have no entry in the
   frequency table). A `*_cut_ga_thres` embedding variant exists in the development repository and
   would be the matched choice for the default annotation settings.
4. **`--embedding_method mean` exports embeddings before NaN filtering**, `sif` after (and after PC
   removal). `bgc_embeddings-mean.tsv` and `bgc_embeddings-sif.tsv` therefore do not have the same row sets.
5. **Unannotated BGCs get string cluster ids** while BIRCH-assigned ones are integers; consumers of
   `bgc_clust-*.tsv` should coerce the column to a single type.
6. **Dereplication is coordinate-based only.** Identical BGCs occurring on *different* contigs
   (e.g. the same locus in two MAGs, or a locus split across two contigs by the assembler) are not
   merged; they are treated as distinct catalog entries and are expected to be reunited later at
   the GCF level.
7. *(Fixed after the Nextflow port, with §11.16.)* **Case 13 in `dereplicate_shared_contigs`**
   (`bin/utils.py`, formerly `dsc.py`) deleted `dereplicated_bgcs[tool2][id1]` (`id1`, not `id2`)
   when the two intervals are identical. It now deletes `id2`.
8. *(Resolved by the Nextflow port.)* The old `dereplicate.py` `--metadata` path crashed with a
   `TypeError` when used without `--input_dir`. `bin/sedna-module-4.py` takes only `--input_dirs`
   and a required `--output_dir`.
9. **`map_synonyms` does an unguarded dict lookup and will crash on an unmapped class.**
   `synonyms_dicts[row_tool][bgc_class_i]` (`bin/utils.py::map_synonyms`) raises `KeyError` whenever a
   tool reports a class absent from its column of `bgc_class_synonyms.tsv`, which aborts the whole
   MODULE4 run at the last step — after the GBKs have already been copied. The risk is
   not hypothetical: commit `ac44bcf` ("added siderophore for antiSMASH — although this was
   supposed to be obsolete it can appear") was exactly this failure. GECCO is the most exposed
   tool, with only 8 label forms in the table against a vocabulary that can emit other
   `;`-joined combinations. A `.get(bgc_class_i, bgc_class_i)` — or an explicit `unknown` —
   would degrade instead of aborting.
10. **The new `bgc_class_formatted` column leaks into the audit tables, unevenly.** `dereplicate_shared_contigs` stores
    *the same* pandas Series object in `dereplicated_bgcs` and in the overlapped / partially
    overlapped / non-overlapped dictionaries, and `row_as_df.loc['bgc_class_formatted'] = …`
    enlarges that Series in place (verified). So every audit row that also won dereplication gains
    the column while the others do not, and `create_df`'s `pd.concat` fills the rest with `NaN`.
    The dereplicated catalog itself is correct; only the audit tables are inconsistent.
11. **`bgc_class_synonyms.tsv` has duplicate keys per tool column.** `NRPS` and `NRPS-like` each
    appear on several rows (differing only in the GECCO/deepBGC columns), so
    `set_index(tool_column)[ref_column].to_dict()` silently keeps the **last** occurrence. Harmless
    today because the duplicates all map to the same `bgc_class_formatted` value, but it means the
    table cannot express a tool-pair-dependent mapping, and a future edit that breaks that
    coincidence would fail silently. Note also that pandas reads the literal `NA` cells as `NaN`,
    which collapses all of them into a single `NaN` key.
12. **The upstream attribution comment on `format_gbks` was removed.** `bin/src/bgc_annot/utilities.py` (now `bin/utils.py` §3.11)
    previously carried "Adapted from Satria A. Kautsar / Wageningen University & Research /
    Copyright (C) 2020"; it was deleted in the same period the repository gained a CC BY-NC-SA 4.0
    `LICENSE`. If that function is still derived from BiG-SLiCE's code, removing the notice does not
    remove the obligation, and the two changes together are worth a look before release — the
    upstream licence may also be incompatible with relicensing under CC BY-NC-SA.
13. *(Fixed in the Nextflow port.)* **The README contained a stray `test push` line** between the beta warning and the workflow
    figure, apparently left over from verifying push access.

14. *(Fixed in the Nextflow port.)* **`map_synonyms` looked up the wrong tool's column after the
    second pairwise merge.** `dereplicate_shared_contigs` takes `tool1 = metadata_df1['tool'].iloc[0]`,
    and from the second merge round on `metadata_df1` is the merged catalog. Its rows (e.g. deepBGC
    BGCs) are therefore filed under the key of the catalog's *first* row's tool (e.g. `antismash`).
    `map_synonyms` selected the synonyms column from that dict key, so a deepBGC `Polyketide` was
    looked up in `bgc_class_antismash` and raised `KeyError`. The legacy pipeline aborts on the shipped
    test sample `OceanDNA-b45222_redu` for this reason. `map_synonyms` now selects the column from each
    row's own `tool` field, as §3.3 describes. The mis-filing itself is left in place, which keeps the
    port output-identical otherwise. It only affects which dict key a row sits under, not the row's
    content or the written tables.
15. *(Fixed after the Nextflow port.)* **antiSMASH and deepBGC coordinates were relative to the BGC,
    not the contig.** `antismash_annot_parser` and `deepbgc_annot_parser` read `start`/`end` from a
    feature inside the per-BGC GBK, which holds only the BGC's sequence, so `start` was ≈ 0. GECCO's
    coordinates (from its TSV) were already on the contig, so module 4 compared positions on different
    frames. The parsers now add the BGC's offset on the contig (antiSMASH `Orig. start`, deepBGC
    accession suffix `_<start>-<end>`). This changes the coordinates and `bgc_id` of every antiSMASH or
    deepBGC BGC not starting at position 0, and can change dereplication results. Details in
    `bugs/001-antismash-deepbgc-relative-coordinates.md`.
16. *(Fixed after the Nextflow port.)* **A BGC that lost a comparison could come back.**
    `dereplicate_shared_contigs` compares every pair of BGCs of a shared contig. The disjoint cases
    6 and 12 add both BGCs to the catalog, so a BGC already superseded by a longer one in an earlier
    pair was re-added whenever a later pair found it disjoint from another BGC. Any contig with two
    or more predictions from one tool was exposed, regardless of the number of tools. Every loser is
    now recorded and removed once all pairs of the round are compared; the 13 cases, longest-wins
    and `--overlap_thres` are unchanged. Details in `bugs/002-dereplication-readds-superseded-bgcs.md`.
---

## 12. References

- **antiSMASH** — Blin et al., *Nucleic Acids Res.* (antiSMASH 6). <https://github.com/antismash/antismash>
- **deepBGC** — Hannigan et al., *Nucleic Acids Res.* 2019, "A deep learning genome-mining strategy
  for biosynthetic gene cluster prediction". <https://github.com/Merck/deepbgc>
- **GECCO** — Carroll, Larralde et al., "Accurate de novo identification of biosynthetic gene
  clusters with GECCO". <https://github.com/zellerlab/GECCO>
- **pyHMMER** — Larralde & Zeller, *Bioinformatics* 2023. <https://github.com/althonos/pyhmmer>
- **cath-resolve-hits** — Lewis et al., *Bioinformatics* 2019 (cath-tools).
  <https://github.com/UCLOrengoGroup/cath-tools>
- **Pfam** — Mistry et al., *Nucleic Acids Res.* 2021.
- **MIBiG 3.1** — Terlouw et al., *Nucleic Acids Res.* 2023.
- **word2vec** — Mikolov et al., 2013. R implementation: the `word2vec` CRAN package.
- **SIF** — Arora, Liang & Ma, ICLR 2017. <https://openreview.net/pdf?id=SyK00v5xx>
- **BIRCH** — Zhang, Ramakrishnan & Livny, *SIGMOD* 1996.
- **BiG-SCAPE / CORASON** — Navarro-Muñoz et al., *Nat. Chem. Biol.* 2020,
  doi:10.1038/s41589-019-0400-9.
- **BiG-SLiCE** — Kautsar et al., *GigaScience* 2021.
- **OceanDNA catalog** (test data) — Nishimura & Yoshizawa, *Sci. Data* 2022.

---

## 13. Nextflow restructuring (implemented)

The pipeline was ported from the `sedna.sh` dispatcher plus hand-written loops to Nextflow,
following the sibling **MG-Clust** pipeline (`~/workspace/repos/tools/MG-Clust`). §1.2 describes
the resulting layout. This section records the decisions behind it.

### 13.1 Why

The previous design chained modules by hand: the README's getting-started section was a sequence of
`for s in ${SAMPLES}; do sedna.sh --module … ; done` loops. That left no dependency graph, no
resumability, no per-process resource declarations and no parallelism beyond what the user wrote.
A workflow manager provides all of these, and MG-Clust already demonstrated the pattern in this
codebase family.

### 13.2 Module mapping

The original design proposal had four modules, with the three predictors inside one. They need three images,
so each predictor became its own module, giving six:

| Process | Former module(s) | Granularity |
|---|---|---|
| `MODULE1` | `run_antismash` | per sample |
| `MODULE2` | `run_deepbgc` | per sample |
| `MODULE3` | `run_gecco` | per sample |
| `MODULE4` | `dereplicate` | per sample |
| `MODULE5` / `MODULE5_PRECOMPUTED` | `annot_cds` / `import_annot_cds` | per sample |
| `MODULE6` | `cluster_sif` / `cluster_mean` (`--embedding_method sif\|mean`) | once, over all samples |
| *(no process)* | `get_bgc_coverage`, `query_ai` | standalone `bin/` tools |

Modules 1–3 run concurrently instead of sequentially as `run_all` did.

### 13.3 The four-layer module pattern

Every module consists of a Python script (`bin/sedna-module-<N>.py`), a Nextflow process
(`modules/sedna-module-<N>.nf`), a Dockerfile (`docker/Dockerfile.module-<N>`) and a conda env
(`docker/resources/requirements-module-<N>.yml`). A new module needs all four, plus an `include`
and a process call in `sedna.nf` and a `docker build` line in `docker/dockerbuild_commands.sh`.

The former `bin/src/bgc_annot` and `bin/src/bgc_derep` packages were merged into `bin/utils.py`,
which every script imports. The module-2 image is Python 3.7, so `utils.py` imports only the stdlib
at module level, imports pandas/biopython inside the functions that use them, and avoids syntax
newer than 3.7.

### 13.4 Points resolved during the port

- **`task.cpus`, not `params.nslots`.** Module `.nf` files pass `--nslots ${task.cpus}`, so
  `nextflow.config` is the single source of truth for resources. The scripts' former `--threads`
  argument is now `--nslots`.
- **Databases stay bind-mounted.** Pfam-A, the antiSMASH databases and the deepBGC models live under
  `params.db_dir` (`~/.sedna/db`) and are mounted through `containerOptions`. `install/get_databases.sh`
  became `bin/download_dbs.py`. Module 1 now passes `--databases` to antiSMASH, which the old wrapper
  never did, and module 2 exports `DEEPBGC_DOWNLOADS_DIR`.
- **Isolated work dirs.** Module 4 copies the catalog GBKs instead of symlinking them, writes `link`
  as a relative path, and remaps each per-tool `file` path onto the staged copy (§3.4).
- **Deterministic merge order.** Module 4 merges antiSMASH → deepBGC → GECCO, where the old script
  used filesystem walk order.
- **Deterministic GCFs.** BIRCH depends on insertion order. The old `cluster_sif.py` used `os.walk`
  order, which is filesystem-dependent, and a naive port would use Nextflow task-completion order, so
  the GCFs could change between runs on identical input. MODULE6 now sorts BGCs by (sample name, file
  name). As a result its cluster ids can differ from a legacy run on the same data. Fed the legacy row
  order, it reproduces the legacy clustering exactly.
- **Empty samples.** A sample with no predicted BGC yields an empty catalog instead of aborting the
  run. Only module 6 fails, and only if no sample had any BGC.
- **Memory.** The `memory` directive is declared per process and scales with `task.attempt` under
  the `slurm` profile. Unlike MG-Clust's MEGAHIT/bbduk/DuckDB, none of SeDNA's tools sizes itself
  from host RAM, so no `--max_mem` threading was needed.
- **Scientific logic unchanged.** The port is structural. The caveats in §11 (SIF weight fallback,
  case 13, unmapped classes in `map_synonyms`, …) are deliberately carried over unchanged. There is one
  exception: §11.14, a crash that stopped the legacy pipeline from completing on its own test data.
- **`sedna.sh` is legacy.** It and the conda installer are kept in `old/` as a reference
  implementation.
