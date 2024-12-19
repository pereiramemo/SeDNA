#!/usr/bin/env python3

###############################################################################
## 1. Set env
###############################################################################

import argparse
import pandas as pd
import numpy as np
import re
import os
import sys
import shutil
import pickle
from collections import Counter
from sklearn.cluster import Birch
from sklearn.decomposition import PCA

current_dir = os.path.dirname(os.path.abspath(__file__))
resources_dir=os.path.join(os.path.dirname(current_dir),"resources")
home_directory = os.getenv('HOME')

###############################################################################
## 2. Parse input data and optional arguments
###############################################################################

parser = argparse.ArgumentParser(prog='cluster.py', \
                                 description='Cluster BGCs into GCFs.')

# general parameters
parser.add_argument("--input_dir", help="Input directory containing the .*_annotdoms_resolved.tsv files.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--threshold", default = 1, help="BIRCH clustering threshold")
parser.add_argument("--embeddings_file", default = f"{home_directory}/.local/share/sedna/embeddings/dom_embeddings_mtx_mibig_pfam.pkl", help="Domain embeddings pkl file.")
parser.add_argument("--counts_file", default = f"{home_directory}/.local/share/sedna/embeddings/domain2counts_mibig3.1_vs_pfam.pkl", help="Domain counts pkl file.")
parser.add_argument("--alpha", default = 0.001, help="SIF a parameter. Default: 0.001.")
parser.add_argument("--bgc_embeddings_tsv", default = None, help="BGC embeddings output table (tsv)")
parser.add_argument("--output_tsv", default = "bgc_clust_output.tsv", help="Output clustering table (tsv).")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")

# Get general parameters
args = parser.parse_args()
input_dir = args.input_dir
sample_name = args.sample_name
threshold = float(args.threshold)
alpha = float(args.alpha)
embeddings_file = args.embeddings_file
counts_file = args.counts_file
bgc_embeddings_tsv = args.bgc_embeddings_tsv
output_tsv = args.output_tsv
overwrite = args.overwrite

###############################################################################
## 3. Check if output_tsv exists
###############################################################################

if os.path.exists(output_tsv) is True:
    
    if args.overwrite is False:
        print(f'Output dir {output_tsv} already exists. Use --overwrite to overwrite')
        sys.exit(0)
    
    if args.overwrite is True:
        try:
            os.remove(output_tsv)
        except OSError as e:
            print(f'Error: {e}')
        
###############################################################################
## 4. Load data from resources
###############################################################################

# Load embeddings
with open(embeddings_file, 'rb') as f:
    domain_embeddings = pickle.load(f)

# Load domain counts
with open(counts_file, 'rb') as f:
    dom2counts = pickle.load(f)

dom2counts_dict = dom2counts.set_index('domain_id')['count'].to_dict()

for dom in list(dom2counts_dict.keys()):
    if dom2counts_dict[dom] < 10:
        del dom2counts_dict[dom]

total_words = sum(dom2counts_dict.values())
weights = {word: alpha / (alpha + (freq / total_words)) for word, freq in dom2counts_dict.items()}

###############################################################################
## 5. Create functions
###############################################################################

def map_domain2embeddings(domains, domain_embeddings, weights):
  
    """
    Compute the SIF-weighted sentence (BGC) embedding for a single sentence.
    
    Args:
        domains (list): A single sentence (BGC).
        domain_embeddings (dict): Pre-trained domain embeddings.
        weights (dict): SIF weights for domains.
    
    Returns:
        dict: A dictionary containing the weighted sentence embedding and a list of weighted domain embeddings.
    """
    
    vector_size = len(next(iter(domain_embeddings.values())))
    embedding_mean = np.zeros(vector_size)
    total_weight = 0

    for domain in domains:
        if domain in domain_embeddings:
          
            weight = weights.get(domain, 1e-3)
            embedding_weighted = weight * np.array(domain_embeddings[domain])
            
            embedding_mean += embedding_weighted    
            
    embedding_norm = np.linalg.norm(embedding_mean)         
    embedding_mean = embedding_mean / embedding_norm 
        
    return(embedding_mean)


def remove_pc(embeddings):
    """
    Remove the first principal component from the sentence embeddings.
    
    Args:
        embeddings (np.ndarray): Array of sentence embeddings.
    
    Returns:
        np.ndarray: Adjusted sentence embeddings with the first principal component removed.
    """
    pca = PCA(n_components=1)
    pca.fit(embeddings)
    pc = pca.components_
    output = embeddings - embeddings.dot(pc.T) * pc
    return output

###############################################################################
## 6. Load input data
###############################################################################

# Initialize list of files.
matching_files = []
for root, dirs, files in os.walk(input_dir):
    for file in files:
        if re.search(".*_annotdoms_resolved.tsv", file):
            if os.stat(os.path.join(root, file)).st_size > 0:
                matching_files.append(os.path.join(root, file))

# Initialize dictionary of lists, echa list containing BGC domains.
bgc_domains = {}
# Load BGCs
for file in matching_files:
    bgc = os.path.basename(file).replace("_annotdoms_resolved.tsv", "")
    df_tmp = pd.read_csv(file, sep=' ', comment='#', header=None)
    bgc_domains[bgc] = df_tmp.iloc[:, 1].tolist()
    
###############################################################################
## 7. Map domains to embeddings and compute weighted mean
###############################################################################

bgc_mean_embeddings = {}

for bgc in bgc_domains:

    # Deduplicate elements in list
    bgc_domains[bgc] = list(dict.fromkeys(bgc_domains[bgc]))     
    bgc_mean_embeddings[bgc] = map_domain2embeddings(bgc_domains[bgc], 
                                                     domain_embeddings, 
                                                     weights)
                                                  
###############################################################################
## 8. Convert list of weighted mean embeddings to df and np objects
###############################################################################

bgc_mean_embedding_df = pd.DataFrame(bgc_mean_embeddings).T
        
###############################################################################
# 9. Identify and remove rows with all NaNs
###############################################################################

i = bgc_mean_embedding_df.isna().all(axis=1)
bgcs_ids_unannot = list(bgc_mean_embedding_df[i].index)
bgc_mean_embedding_df = bgc_mean_embedding_df[~i]
        
###############################################################################
## 10. Remove first PC 
###############################################################################

# convert to np to run PCA in remove_pc
bgc_mean_embedding_np = bgc_mean_embedding_df.to_numpy()
# remove fist PC
bgc_mean_embedding_minus_pc1_np = remove_pc(bgc_mean_embedding_np)
# convert back to df
bgc_mean_embedding_minus_pc1_df = pd.DataFrame(bgc_mean_embedding_minus_pc1_np,
                                               index=bgc_mean_embedding_df.index,
                                               columns = bgc_mean_embedding_df.columns
                                              )

if bgc_embeddings_tsv is not None:
    bgc_mean_embedding_minus_pc1_df.to_csv(bgc_embeddings_tsv, sep='\t', index=True)

###############################################################################
## 11. Cluster BGCs
###############################################################################

birch = Birch(
            n_clusters=None,  
            compute_labels=False, 
            copy=False 
        )
        
birch.threshold = threshold

birch.branching_factor = bgc_mean_embedding_minus_pc1_df.shape[0]
birch.fit(bgc_mean_embedding_minus_pc1_df)
brc_cluster = birch.predict(bgc_mean_embedding_minus_pc1_df)

###############################################################################
## 12. Format clustering as df
###############################################################################

bgc_ids = bgc_mean_embedding_minus_pc1_df.index

output_df = pd.DataFrame({
    'bgc_id': bgc_ids,
    'cluster_id': brc_cluster
})

###############################################################################
## 13. Add un annot BGCs
###############################################################################

cluster_id_max = output_df['cluster_id'].max()
cluster_ids_unannot = list(range(cluster_id_max +1, 
                                 cluster_id_max + len(bgcs_ids_unannot) +1))

output_unannot_df = pd.DataFrame({'bgc_id': bgcs_ids_unannot, 
                                  'cluster_id': cluster_ids_unannot})

output_df = pd.concat([output_df, output_unannot_df], ignore_index=True)

###############################################################################
## 14. Export clustering
###############################################################################

output_df.to_csv(output_tsv, sep='\t', index=False)

print("Clustering executed successfully")
