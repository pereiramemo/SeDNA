# SeDNA Bioprospecting Pipeline
Sedna is a bioinformatic pipeline dedicated to the bioprospecting analysis of Biosynthetic Gene Clusters in metagenomic data.

The tasks performed by Sedna are the following.  
**Identify BGC sequences**. To take full advantage of current and future developments in BGC prediction, our pipeline implements a plug-and-play integration of third-party BGC annotation tools. This allows us to use an ensemble approach to improve both sensitivity and specificity in BGC detection. At present, the tools included are antiSMASH, gecco, and deepBGC.  
**Creating a BGC catalog**. Since we integrate different BGC prediction tools, some BGC sequences will be predicted by more than one tool, leading to duplications in our results. To generate a catalog of unique BGC sequences, our pipeline implements a dereplication step. This process uses the contig IDs and coordinates of the BGC sequences to remove duplicates. Accordingly, this approach enables us to obtain a catalog of unique BGC sequences that integrates predictions from multiple BGC annotation tools.  

# How to install

SeDNA is a [Nextflow](https://www.nextflow.io/) pipeline. Each module runs in its own container, so the only host requirements are:

- Nextflow (≥ 24.x),
- Docker (default `standard` profile), or Apptainer (`slurm` profile),
- Python 3 (for the two setup scripts in `bin/`).

1. Clone the repository:
```
git clone https://github.com/new-atlantis-labs/SeDNA.git
cd SeDNA
```

2. Build the four module images. This must be run from the repository root:
```
bash docker/dockerbuild_commands.sh
```
On a cluster with Apptainer, cache the images on the login node instead, once they are available in the registry:
```
bin/download_imgs.py
```

3. Download the databases (antiSMASH, deepBGC) into `~/.sedna/db`:
```
bin/download_dbs.py
```
If the databases were already fetched by the old `install/get_databases.sh`, you can reuse them with `--db_dir ~/.local/share` instead.

# Getting started

Print the full list of parameters:
```
nextflow run sedna.nf --help
```

SeDNA reads a tab-separated samplesheet with one row per assembly or MAG (FASTA):

```
sample_name	fasta
OceanDNA-b11979_redu	tests/data/oceandna/OceanDNA-b11979_redu.fasta
OceanDNA-b45222_redu	tests/data/oceandna/OceanDNA-b45222_redu.fasta
OceanDNA-b9899_redu	tests/data/oceandna/OceanDNA-b9899_redu.fasta
```

The file above ships as `tests/data/from_fasta.tsv`. It lists three toy MAGs from the [OceanDNA catalog](https://www.nature.com/articles/s41597-022-01392-5), artificially subsampled to include only a few BGCs. To run the whole pipeline on them:

```
nextflow run sedna.nf \
    --from_fasta_tsv tests/data/from_fasta.tsv \
    --output_dir tests/sedna-output
```

Add `-resume` to continue an interrupted run, or to re-run with different parameters: only the modules affected by a changed parameter, and those downstream of them, are re-run. For example, changing `--overlap_thres` re-runs only module 4. Add `-profile slurm` to submit to a SLURM cluster with Apptainer.

The pipeline runs four modules:

| Module | Task | Runs |
|---|---|---|
| 1 | BGC annotation with [antiSMASH](https://github.com/antismash/antismash) | per sample |
| 2 | BGC annotation with [deepBGC](https://github.com/Merck/deepbgc) | per sample |
| 3 | BGC annotation with [GECCO](https://github.com/zellerlab/GECCO) | per sample |
| 4 | Dereplication of the three tools' predictions | per sample |

Modules 1–3 run concurrently. Any of them can be dropped with `--skip_antismash`, `--skip_deepbgc` or `--skip_gecco`. Since the same sequences are annotated with different tools, some BGCs are predicted more than once; module 4 merges these duplicates into one catalog per sample, using the contig IDs and coordinates of the predictions. It works with any number of annotation tools, and writes the catalog in three versions (`longest`, `core`, `extended`; see Outputs). `--stop_at_module N` stops the pipeline after module N.

### Outputs

```
tests/sedna-output/
  intermediate/
    module-1/<sample>/          antiSMASH: gbks/, annot_metadata.tsv
    module-2/<sample>/          deepBGC:   gbks/, gbks_removed/, annot_metadata.tsv
    module-3/<sample>/          GECCO:     gbks/, annot_metadata.tsv
    module-4/<sample>/          dereplicated/{longest,core,extended}/ (the BGC catalog),
                                overlapped/, partially_overlapped/, non_overlapped/
```

The de-replicated catalog comes in three versions of the same BGCs, which differ only in the coordinates of BGCs predicted by several tools:

| Version | Coordinates |
|---|---|
| `longest` | the longest prediction |
| `core` | the region supported by the most tools (their intersection when they all overlap) |
| `extended` | from the leftmost start to the rightmost end of all predictions |

Each `module-4/<sample>/dereplicated/<version>/annot_metadata.tsv` has the fields `acc`, `bgc_class`, `start`, `end`, `length`, `on_edge`, `contig_id`, `contig_length`, `file`, `sample_name`, `tool`, `bgc_id`, `link`, `bgc_class_formatted`, `tools`, `n_tools` and `members` (plus `core_support` in `core`). The BGC sequences themselves are in `dereplicated/<version>/gbks/` as GenBank files named `<bgc_id>.gbk`. `bgc_id` is `<tool>__<sample>__<contig>__<start>__<end>`, where `<tool>` lists all supporting tools joined by `+` in `core` and `extended`.

Outputs of modules 1–3 are published only with `--full_output true` (the default).

# How to uninstall

Remove the images (`docker rmi ghcr.io/pereiramemo/sedna/module-{1..4}:latest`) and the database and cache directory `~/.sedna`.


# Citation

If you use Sedna in your research, please cite:

Pereira Flores, E. (2026). *pereiramemo/SeDNA: v1.0.0* (beta). Zenodo. [https://doi.org/10.5281/zenodo.23224099](https://doi.org/10.5281/zenodo.23224099)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23224099.svg)](https://doi.org/10.5281/zenodo.23224099)

```bibtex
@software{pereira_flores_2026_sedna,
  author    = {Pereira Flores, Emiliano},
  title     = {pereiramemo/SeDNA: v1.0.0},
  year      = {2026},
  publisher = {Zenodo},
  version   = {beta},
  doi       = {10.5281/zenodo.23224099},
  url       = {https://doi.org/10.5281/zenodo.23224099}
}
```

# License

Sedna is licensed under the [Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License](https://creativecommons.org/licenses/by-nc-sa/4.0/) (CC BY-NC-SA 4.0). See the [LICENSE](LICENSE) file for the full text.

Sedna was developed at NewAtlantis Labs, which has since been acquired by Ocean BioMetrics.

For licensing questions, please contact Emiliano Pereira at [pereiramemo@gmail.com](mailto:pereiramemo@gmail.com).
