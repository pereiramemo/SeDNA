# Bioprospecting-Internal
BGCminer is a bioinformatic pipeline dedicated to the bioprospecting analysis of Biosynthetic Gene Clusters in metagenomic data (see fig. 1 below).

![Figure 1: BGCminer workflow](https://github.com/new-atlantis-labs/Bioprospecting-Internal/blob/main/images/bgc_miner.png)

The tasks performed by BGCminer are the following.  
**Identify BGC sequences**. To take full advantage of current advancements and future developments in BGC prediction, our pipeline implements a plug-and-play integration of third-party BGC annotation tools. This allows us to use an ensemble approach to improve both sensitivity and specificity in BGC detection. At present, the tools included antiSMASH, gecco.  
**Creating a BGC catalog**. Since we integrate different BGC prediction tools, some BGC sequences will be predicted by more than one tool, leading to duplications in our results. To generate a catalog of unique BGC sequences, our pipeline implements a dereplication step. This process uses the contig IDs and coordinates of the BGC sequences to remove duplicates. Accordingly, this approach enables us to obtain a catalog of unique BGC sequences that integrates predictions from multiple BGC annotation tools.  
**BGC clustering**. To organize the large diversity of BGC sequences found in microbial communities, we implemented a novel clustering methodology to group BGC sequences responsible for producing similar compounds into Gene Cluster Families (GCFs).  
**Functional annotation**. By integrating BGC reference databases, such as MIBiG, into our clustering methodology, the pipeline provides a functional and novelty assessment of environmental BGC sequences. Additionally, the GCF composition can be used to estimate the biosynthetic diversity in metagenomic samples. These analyses help identify candidate BGC sequences for specific industrial applications and facilitate the selection of taxa or environments with high biosynthetic novelty and diversity.  
