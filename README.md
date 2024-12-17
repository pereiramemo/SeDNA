# SeDNA Bioprospecting Pipline
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
SeDNA is composed of different modules that must be run sequentially in a specific order.  
First, we need to annotate the BGC sequences in the assembled metagenomics sequences of Metagenome-Assembled-Genomes (all as fasta files). In this example, we will be analyzing three toy MAGs from the [OceanDNA catalog](https://www.nature.com/articles/s41597-022-01392-5), which can be found in the `tests/data` folder with the repository. 
The current version of SeDNA integrates three BGC annotation tools: [antiSMASH](https://github.com/antismash/antismash), [deepBGC](https://github.com/Merck/deepbgc), and [GECCO](https://github.com/zellerlab/GECCO). We can run these one-by-one, using the `run_<tool>` modules or all at once with the `run_all` module.
Let’s navigate to the root of the repository, and get all the samples in a variable:

```
SAMPLES=$(ls "tests/data/OceanDNA"-*.fa)
```

Now, we iterate over these samples, and annotate the BGC sequences utilizing the three tools. For simplicity we will use the `run_all` module. 
```
for s in ${SAMPLES}; do

SAMPLE_NAME=$(basename "${s}" .fa)
    echo "${SAMPLE_NAME}"
    sedna.sh --module run_all \
        --params "--input_sample ${s} \
                  --sample_name ${SAMPLE_NAME} \
                  --output_dir ${REPO_DIR}/test/output/sedna_output";

done
```











