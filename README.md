# SeDNA Bioprospecting Pipline
Sedna is a bioinformatic pipeline dedicated to the bioprospecting analysis of Biosynthetic Gene Clusters in metagenomic data (see fig. 1 below).

![Figure 1: sedna workflow](https://github.com/new-atlantis-labs/Bioprospecting-Internal/blob/main/images/sedna.png)

The tasks performed by Sedna are the following.  
**Identify BGC sequences**. To take full advantage of current and future developments in BGC prediction, our pipeline implements a plug-and-play integration of third-party BGC annotation tools. This allows us to use an ensemble approach to improve both sensitivity and specificity in BGC detection. At present, the tools included antiSMASH, gecco, and deepBGC.  
**Creating a BGC catalog**. Since we integrate different BGC prediction tools, some BGC sequences will be predicted by more than one tool, leading to duplications in our results. To generate a catalog of unique BGC sequences, our pipeline implements a dereplication step. This process uses the contig IDs and coordinates of the BGC sequences to remove duplicates. Accordingly, this approach enables us to obtain a catalog of unique BGC sequences that integrates predictions from multiple BGC annotation tools.  
**BGC clustering**. To organize the large diversity of BGC sequences found in microbial communities, we implemented a novel clustering methodology to group BGC sequences responsible for producing similar compounds into Gene Cluster Families (GCFs).  
**Functional annotation**. By integrating BGC reference databases, such as MIBiG, into our clustering methodology, the pipeline provides a functional and novelty assessment of environmental BGC sequences. Additionally, the GCF composition can be used to estimate the biosynthetic diversity in metagenomic samples. These analyses help identify candidate BGC sequences for specific industrial applications and facilitate the selection of taxa or environments with high biosynthetic novelty and diversity.  


# How to install

### Recommended pre-installation steps

1. Make sure conda and mamba are installed. 
```
conda --version
```
```
mamba --version
```

If conda is not installed, you can install if following the instrctions [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
If manba is not insalled, you can install it with the following command (after installing conda).
```
conda install conda-forge::mamba
```

2. Clean and update conda
```
conda clean --all -y 
```
```
conda update -n base --all -y
```

### Installation steps

1. Clone the repository  
```
git clone https://github.com/new-atlantis-labs/SeDNA.git
```

2. Navigate to the repositry  
```
cd SeDNA
```

3. Run the `install.sh` script within the `install` direcotry  
```
./install/install.sh
```

4. Run the `get_databases.sh` script within the `install` direcotry  
```
./install/get_databases.sh
```
