# SeDNA — Technical Documentation

**Version documented:** 0.1.0 (`bin/VERSION.txt`)
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
FASTA (MAGs / assembled contigs)
        |
        |  [1] BGC annotation          modules: run_antismash | run_deepbgc | run_gecco | run_all
        |      antiSMASH + deepBGC + GECCO, run independently on the same input
        v
   per-tool GBKs + annot_metadata.tsv  (bgc_annot/inter/<tool>/)
        |
        |  [2] Dereplication           module: dereplicate
        |      coordinate-based merging of the three per-tool predictions
        v
   non-redundant BGC catalog          (bgc_annot/sorted/dereplicated/)
        |
        |  [3] Domain annotation       module: annot_cds  (or import_annot_cds)
        |      pyHMMER hmmsearch vs Pfam-A  ->  cath-resolve-hits
        v
   *_annotdoms_resolved.tsv  (one domain architecture string per BGC)
        |
        |  [4] Embedding + [5] clustering   module: cluster_sif (or cluster_mean)
        |      word2vec domain vectors -> SIF-weighted BGC vector -> BIRCH
        v
   bgc_embeddings.tsv  +  clust.tsv  (bgc_id -> GCF id)
```

### 1.2 Architecture and execution model

`bin/sedna.sh` is a thin dispatcher. Each module is a standalone Python script in `bin/` that is
executed inside its **own conda environment**, activated on the fly by the dispatcher:

| Module | Script | Environment |
|---|---|---|
| `run_antismash` | `bin/run_antismash.py` | `sedna_antismash_env` |
| `run_deepbgc` | `bin/run_deepbgc.py` | `sedna_deepbgc_env` |
| `run_gecco` | `bin/run_gecco.py` | `sedna_gecco_env` |
| `run_all` | the three above, in sequence | (all three) |
| `dereplicate` | `bin/dereplicate.py` | `sedna_general_purpose_env` |
| `annot_cds` | `bin/annot_cds.py` | `sedna_annot_cds_env` |
| `import_annot_cds` | `bin/import_annot_cds.py` | `sedna_import_annot_cds_env` |
| `get_bgc_coverage` | `bin/get_bgc_coverage.py` | `sedna_get_bgc_coverage_env` |
| `cluster_mean` | `bin/cluster_mean.py` | `sedna_general_purpose_env` |
| `cluster_sif` | `bin/cluster_sif.py` | `sedna_general_purpose_env` |

This per-module isolation is what makes the "plug-and-play" integration of third-party BGC
predictors practical: antiSMASH 6.1.1 (Python 3.9), deepBGC 0.1.31 (Python 3.7 + TensorFlow 1.15)
and GECCO 0.9.10 (Python 3.12) have mutually incompatible dependency trees and simply cannot
coexist in a single environment.

Usage:

```bash
conda activate sedna_main_env
sedna.sh --module <module> --params "<module arguments>"
sedna.sh --module <module> --params "--help"     # per-module help
```

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
`--output_dir`, `--threads`, `--overwrite`) and validate the input FASTA before running.

**antiSMASH** (`bin/run_antismash.py`)

```
antismash --cpus <threads> --genefinding-tool prodigal-m --taxon bacteria \
          --allow-long-headers --minlength 1000 --minimal \
          --output-dir <out> <input.fasta>
```

- `prodigal-m` (Prodigal metagenomic/anonymous mode) is the default gene finder — appropriate for
  fragmented metagenomic contigs of unknown provenance.
- `--minimal` skips the optional comparative/analysis modules; only detection is needed here.
- Outputs are the `*.region*.gbk` files, moved into `bgc_annot/inter/antismash/gbks/`.

**deepBGC** (`bin/run_deepbgc.py`)

```
deepbgc pipeline --detector deepbgc --classifier product_class \
                 --prodigal-meta-mode --min-nucl 1000 --output <out> <input.fasta>
```

- deepBGC emits a single multi-record GBK, which is split into one file per BGC
  (`utilities.gbk_splitter`) and renamed to the antiSMASH-style `<contig>.region<NNN>.gbk`
  convention so that all tools produce a uniform file layout.
- Predictions are then filtered **after** the run (`utilities.gbk_filter`) using
  `--score_thres` (default **0.80**) and `--cds_count_thres` (default **3** CDS).
  Filtering is deliberately done post-hoc rather than via tool arguments, so that relaxing the
  stringency later does not require re-running the tool; discarded clusters are retained under
  `gbks_removed/`.
- The reported BGC class is the `product_class` with the highest score.

**GECCO** (`bin/run_gecco.py`)

```
gecco run --genome <input.fasta> --cds 3 --threshold 0.80 --jobs <threads> --output-dir <out>
```

- GECCO's `*_cluster_<n>.gbk` outputs are renamed to `.region<NNN>.gbk`, again for uniformity.

### 2.3 Per-tool metadata

Each wrapper ends by parsing its own GenBank output into a common table,
`bgc_annot/inter/<tool>/annot_metadata.tsv`, via a tool-specific parser in
`bin/src/bgc_annot/utilities.py` (`antismash_annot_parser`, `deepbgc_annot_parser`,
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
- BGC class vocabularies differ between tools; `resources/bgc_class_synonyms.tsv` (75 entries)
  maps antiSMASH / GECCO / deepBGC class labels onto a single harmonised
  `bgc_class_formatted` vocabulary.
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

`bin/dereplicate.py` loads every `annot_metadata.tsv` found under `--input_dir` and merges them
**recursively and pairwise**: tables 1 and 2 are merged, the merged catalog becomes the new
table 1, which is then merged with table 3, and so on (`recursive_dereplication`). For each pair:

**(a) Shared contigs** (`bin/src/bgc_derep/dsc.py::dereplicate_shared_contigs`).
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

**(b) Non-shared contigs** (`bin/src/bgc_derep/dnsc.py`). BGCs on contigs seen by only one of the
two tables cannot conflict with anything, so they are carried over directly (case 14).

**(c) Single input table** (`bin/src/bgc_derep/sst.py`). If only one tool was run, every BGC is
passed straight through (case 14) so that downstream modules receive the same directory layout.

### 3.3 Output

`bin/src/bgc_derep/outputs.py` writes four sibling directories under `bgc_annot/sorted/`:

```
sorted/
  dereplicated/           <- the catalog: one entry per unique locus
    annot_metadata.tsv
    gbks/                 <- symlinks to the winning per-tool GBK files
  overlapped/             <- audit: predictions that matched another prediction
  partially_overlapped/   <- audit: staggered matches below --overlap_thres
  non_overlapped/         <- audit: predictions with no counterpart
```

Each retained BGC receives a stable, self-describing identifier:

```
bgc_id = <tool>__<sample_name>__<contig_id>__<start>__<end>
```

and the metadata gains two columns: `bgc_id` and `link` (absolute path of the symlinked GBK).
Because the catalog stores symlinks rather than copies, the provenance of every sequence back to
the tool that produced it is preserved.

---

## 4. BGC domain annotation

Clustering operates on **protein-domain composition**, not on nucleotide or protein sequence
similarity. Each BGC is reduced to the ordered list of Pfam domains encoded in its CDSs — the
"sentence of domains" that the embedding step consumes.

### 4.1 `annot_cds` — pyHMMER + cath-resolve-hits

`bin/annot_cds.py` performs three steps for every `*.gbk` in the dereplicated catalog:

**Step 1 — CDS extraction.** `utilities.extract_cds` writes the translated CDSs of each BGC to a
FASTA file, using `--cds_id` (default `locus_tag`; `protein_id` for MIBiG-style records) as the
sequence identifier.

**Step 2 — HMM search with [pyHMMER](https://github.com/althonos/pyhmmer)** (0.10.x):

```python
pyhmmer.hmmer.hmmsearch(hmms, seqs, domE=evalue_thres, bit_cutoffs="gathering", cpus=threads)
```

- Database: **Pfam-A** (`~/.local/share/pfam/Pfam-A.hmm`, downloaded by `install/get_databases.sh`).
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

`bin/import_annot_cds.py` produces the same `*_annotdoms_resolved.tsv` files from annotation that
already exists elsewhere: given `--annot_metadata` (BGC coordinates), a `--gff` of gene features
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
(installed to `~/.local/share/sedna/embeddings/`), a plain `dict` mapping a Pfam family name to a
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

`bin/cluster_sif.py` implements **Smooth Inverse Frequency** weighting (Arora, Liang & Ma, ICLR
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

The result is written to `--bgc_embeddings_tsv` (rows = BGCs, 100 columns) and is directly reusable
for downstream novelty and diversity analyses.

### 5.4 `cluster_mean` — the simpler alternative

`bin/cluster_mean.py` is the earlier/baseline variant: it maps domains to vectors and takes the
**unweighted arithmetic mean**, with no SIF weighting, no L2 normalisation and no principal
component removal, then clusters with BIRCH exactly as `cluster_sif` does. It is retained for
comparison; `cluster_sif` is the recommended module.

---

## 6. BGC clustering into GCFs — BIRCH

### 6.1 Configuration

Both clustering modules use scikit-learn's
[`Birch`](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.Birch.html)
(Balanced Iterative Reducing and Clustering using Hierarchies, Zhang et al., SIGMOD 1996):

```python
birch = Birch(n_clusters=None, compute_labels=False, copy=False)
birch.threshold = <--threshold>
birch.branching_factor = <number of BGCs>
birch.fit(X)
labels = birch.predict(X)
```

- `n_clusters=None` — **no global clustering step**. The leaf subclusters of the CF-tree *are* the
  GCFs. There is no need to state in advance how many gene cluster families a metagenome contains
  (which is the actual unknown), and no arbitrary cut of a dendrogram.
- `threshold` (`--threshold`, default 1; 0.5 in the README example) is the **only** parameter that
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

`--output_tsv` (default `bgc_clust_output.tsv`) is a two-column table:

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
> re-weighting is currently baked into it, and neither `cluster_sif.py` nor `cluster_mean.py`
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

## 9. Auxiliary module: `get_bgc_coverage`

`bin/get_bgc_coverage.py` estimates read coverage (i.e. relative abundance) for each BGC in the
catalog, from either a BAM file of reads mapped to the assembly or a per-contig coverage table
(e.g. produced by VEBA), using `bedtools`/`pysam` over the BGC coordinates in `annot_metadata.tsv`.
This turns presence/absence GCF profiles into quantitative ones for cross-sample comparison.

---

## 10. Installation and data dependencies

```bash
./install/install.sh        # creates sedna_main_env + one env per module, installs bin/ into each
./install/get_databases.sh  # downloads databases and installs the embedding models
./install/uninstall.sh      # removes all of the above
```

`get_databases.sh` populates:

| Path | Content |
|---|---|
| `~/.local/share/pfam/Pfam-A.hmm` | Pfam-A HMM library (EBI, current release) |
| `~/.local/share/antismash/` | antiSMASH databases (`download-antismash-databases`) |
| (deepBGC default location) | deepBGC models (`deepbgc download`) |
| `~/.local/share/sedna/embeddings/` | `dom_embeddings_mtx_mibig_pfam.pkl`, `domain2counts_mibig3.1_vs_pfam.pkl` |

A `docker/Dockerfile` is also provided.

---

## 11. Implementation notes and caveats

These are behaviours of the current code worth being aware of when interpreting results or
modifying the pipeline.

1. **SIF weight fallback inverts the SIF intent for rare domains.** In
   `cluster_sif.py::map_domain2embeddings`, a domain absent from the filtered weight table receives
   `weights.get(domain, 1e-3)`. The weight table only keeps domains with count ≥ 10 (656 of 3,397),
   so of the 1,262 domains in the embedding vocabulary, **580 rare domains (plus 26 keys with no
   count entry, including word2vec's `</s>` token) receive the 1e-3 fallback** — *lower* than the
   weight of the most frequent domain in the corpus (0.0088). Under SIF, rare domains should
   receive the *highest* weights (up to 0.894 here). The practical effect is that rare — and
   therefore potentially most diagnostic — domains are nearly excluded from the BGC vector.
2. **Sum vs. weighted mean.** The docstring says "weighted mean", but the implementation sums the
   weighted domain vectors and divides by the L2 norm; the accumulated `total_weight` variable is
   unused. This is consistent with the SIF paper up to a scale factor, and the subsequent
   L2 normalisation makes the difference immaterial — but the code and docstring disagree.
3. **Embedding corpus vs. default annotation settings.** The shipped embedding matrix and the
   frequency table were derived from the MIBiG 3.1 Pfam annotation produced with an **E-value**
   threshold, while `annot_cds.py` defaults to **gathering thresholds** (`--cut_ga`). The two
   regimes produce slightly different domain vocabularies (26 embedded domains have no entry in the
   frequency table). A `*_cut_ga_thres` embedding variant exists in the development repository and
   would be the matched choice for the default annotation settings.
4. **`cluster_mean` exports embeddings before NaN filtering**, `cluster_sif` after (and after PC
   removal). The two `--bgc_embeddings_tsv` outputs therefore do not have the same row sets.
5. **Unannotated BGCs get string cluster ids** while BIRCH-assigned ones are integers; consumers of
   `clust.tsv` should coerce the column to a single type.
6. **Dereplication is coordinate-based only.** Identical BGCs occurring on *different* contigs
   (e.g. the same locus in two MAGs, or a locus split across two contigs by the assembler) are not
   merged; they are treated as distinct catalog entries and are expected to be reunited later at
   the GCF level.
7. **Case 13 in `dsc.py`** deletes `dereplicated_bgcs[tool2][id1]` (note: `id1`, not `id2`) when the
   two intervals are identical. Worth reviewing.
8. **`dereplicate.py`'s `--metadata` path is not fully usable on its own.** With `--metadata` but
   no `--input_dir`, the default output directory is still derived from `input_dir`
   (`os.path.abspath(None)`) and raises a `TypeError`; passing neither argument raises a
   `TypeError` from `os.path.exists(None)` instead of the intended error message. In practice
   `--input_dir` (optionally together with `--output_dir`) is the supported invocation.

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
