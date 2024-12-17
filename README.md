# SeDNA Bioprospecting Pipeline
Sedna is a bioinformatic pipeline dedicated to the bioprospecting analysis of Biosynthetic Gene Clusters in metagenomic data (see fig. 1 below).

![Figure 1: sedna workflow](https://github.com/new-atlantis-labs/Bioprospecting-Internal/blob/main/images/sedna.png)

The tasks performed by Sedna are the following.  
**Identify BGC sequences**. To take full advantage of current and future developments in BGC prediction, our pipeline implements a plug-and-play integration of third-party BGC annotation tools. This allows us to use an ensemble approach to improve both sensitivity and specificity in BGC detection. At present, the tools included antiSMASH, gecco, and deepBGC.  
**Creating a BGC catalog**. Since we integrate different BGC prediction tools, some BGC sequences will be predicted by more than one tool, leading to duplications in our results. To generate a catalog of unique BGC sequences, our pipeline implements a dereplication step. This process uses the contig IDs and coordinates of the BGC sequences to remove duplicates. Accordingly, this approach enables us to obtain a catalog of unique BGC sequences that integrates predictions from multiple BGC annotation tools.  
**BGC clustering**. To organize the large diversity of BGC sequences found in microbial communities, we implemented a novel clustering methodology to group BGC sequences responsible for producing similar compounds into Gene Cluster Families (GCFs).  
**Functional annotation**. By integrating BGC reference databases, such as MIBiG, into our clustering methodology, the pipeline provides a functional and novelty assessment of environmental BGC sequences. Additionally, the GCF composition can be used to estimate the biosynthetic diversity in metagenomic samples. These analyses help identify candidate BGC sequences for specific industrial applications and facilitate the selection of taxa or environments with high biosynthetic novelty and diversity.  

# How to install

### Recommended Pre-installation Steps

1. Ensure Conda and Mamba are installed.
Check their versions with the following commands:
```
conda --version
```
```
mamba --version
```


If Conda is not installed, follow the instructions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).  
If Mamba is not installed, you can install it using Conda:
```
conda install conda-forge::mamba
```

2. Clean and update Conda:
```
conda clean --all -y 
```
```
conda update -n base --all -y
```

### Installation Steps

1. Clone the repository:  
```
git clone https://github.com/new-atlantis-labs/SeDNA.git
```

2. Navigate to the repository: 
```
cd SeDNA
```

3. Run the installation script located in the `install` directory:
```
./install/install.sh
```

4. Download required databases using the provided script: 
```
./install/get_databases.sh
```

# Getting started

After installing SeDNA, you first need to activate the main environment.
```
conda activate sedna_main_env
```
Now, you should be able to run the SeDNA wrapping script:
```
sedna.sh --help
```
We get the following help instructions:
```
Usage: ./bin/sedna.sh [-m <module>] [-o <options>] [-v|--version] [-h|--help]
 
Options:
  -m, --module    Specify the module. Available modules: run_antismash run_deepbgc run_gecco run_all dereplicate annot_cds cluster_mean cluster_sif
  -p, --params    Specify parameters to give to each module
  -v, --version   Display the version information
  -h, --help      Display this help message
```

To get the help of an specific module, we can run the following command:
```
sedna.sh --module <module_name> --params "--help"
```
For example, 
```
sedna.sh --module run_antismash --params "--help"
```
And we get the help instructions:
```
usage: run_antismash.py [-h] --input_sample INPUT_SAMPLE [--threads THREADS] [--sample_name SAMPLE_NAME] [--output_dir OUTPUT_DIR] [--overwrite] [--taxon TAXON]
                        [--genefinding_tool GENEFINDING_TOOL] [--minlength MINLENGTH]

Annotates BGC sequences utilizing the antiSMASH tool.

options:
  -h, --help            show this help message and exit
  --input_sample INPUT_SAMPLE
                        Input fasta file.
  --threads THREADS     Number of threads.
  --sample_name SAMPLE_NAME
                        Sample name.
  --output_dir OUTPUT_DIR
                        Output directory.
  --overwrite           Overwrite output directory.
  --taxon TAXON         antiSMASH - {bacteria,fungi} Taxonomic classification of input sequence.
  --genefinding_tool GENEFINDING_TOOL
                        antiSMASH {glimmerhmm,prodigal,prodigal-m,none,error} Specify algorithm used for gene finding: GlimmerHMM, Prodigal, Prodigal Metagenomic/Anonymous mode, or none. The
                        'error' option will raise an error if genefinding is attempted. The 'none' option will not run genefinding.
  --minlength MINLENGTH
                        antiSMASH - Only process sequences larger than <minlength>.
                        
```
SeDNA is composed of different modules that must be run sequentially in a specific order.

First, we need to annotate the BGC sequences in the assembled metagenomic sequences of Metagenome-Assembled Genomes (MAGs), provided as FASTA files. In this example, we will analyze three toy MAGs from the [OceanDNA catalog](https://www.nature.com/articles/s41597-022-01392-5), which can be found in the tests/data folder within the repository.

The current version of SeDNA integrates three BGC annotation tools: [antiSMASH](https://github.com/antismash/antismash), [deepBGC](https://github.com/Merck/deepbgc), and [GECCO](https://github.com/zellerlab/GECCO). These tools can be run individually using the `run_<tool>` modules or all at once using the `run_all` option.

Let’s navigate to the root of the repository and store all the sample files in a variable:

```
SAMPLES=$(ls "tests/data/OceanDNA"-*_redu.fasta)
```
Let's use the `tests` folder to save our example analysis.
```
OUTPUT_DIR=./tests/output/sedna_output
```
Now, we iterate over the input samples and annotate the BGC sequences using the three tools. For simplicity, we will use the `run_all` option.

```
for s in ${SAMPLES}; do

    SAMPLE_NAME=$(basename "${s}" .fasta)
        echo "${SAMPLE_NAME}"
        sedna.sh --module run_all \
            --params "--input_sample ${s} \
                      --sample_name ${SAMPLE_NAME} \
                      --output_dir ${OUTPUT_DIR}";

done
```

"Since we are annotating the same sequences with different tools, it is expected that some BGC sequences will be predicted by more than one tool. This means we may have duplicate BGC predictions. To obtain a de-replicated catalog, we need to run the `dereplicate` module as follows:

```
for s in ${SAMPLES}; do

    SAMPLE_NAME=$(basename "${s}" .fa)
    echo "${SAMPLE_NAME}"
    sedna.sh --module dereplicate \
    --params "--input_dir ${OUTPUT_DIR}/${SAMPLE_NAME}/bgc_annot/inter"
           
done                     
```

"This script will generate a `sorted` folder within each output directory, where we can find the dereplicated catalog.
Let’s take a look at one of these:

```
ls ${OUTPUT_DIR}/OceanDNA-b11979_redu/bgc_annot/sorted/dereplicated/
```

Once we obtain the de-replicated catalog, we will cluster the BGC sequences predicted in the three samples.
To do this, we must annotate the Pfam domains using the `cds_annot` module:

```
for s in ${SAMPLES}; do

      SAMPLE_NAME=$(basename "${s}" .fa)
      echo "${SAMPLE_NAME}"
           sedna.sh --module annot_cds \
           --params "--input_dir ${OUTPUT_DIR}/${SAMPLE_NAME}/bgc_annot/sorted/dereplicated"
done 

```

Finally, based on the domain annotations, we will cluster the BGC sequences using the cluster_sif module, which implements the [Balanced Iterative Reducing and Clustering using Hierarchies (BIRCH)](https://en.wikipedia.org/wiki/BIRCH) and utilizes the [Smooth Inverse Function (SIF)](https://openreview.net/pdf?id=SyK00v5xx) to construct the BGC embeddings.

```
sedna.sh --module cluster_sif \
    --params "--input_dir  ${OUTPUT_DIR} \
              --threshold 3 \
              --bgc_embeddings_tsv ${OUTPUT_DIR}/bgc_embeddings.tsv \
              --output_tsv ${OUTPUT_DIR}/clust.tsv"
              
```              

The file `${OUTPUT_DIR}/bgc_embeddings.tsv` contains the emebdding of each BGC, which can be useful for downstream analysis, such as novelty or diversity assessment. 

The file `${OUTPUT_DIR}/clust.tsv`contains the actual clustering 



