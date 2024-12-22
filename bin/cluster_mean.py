#!/usr/bin/env python3

###############################################################################
## 1. Set env
###############################################################################

import argparse
import pandas as pd
import re
import os
import sys
import shutil
import pickle
from sklearn.cluster import Birch

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
parser.add_argument("--output_tsv", default = "bgc_clust_output.tsv", help="Output clustering table (tsv).")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")

# Get general parameters
args = parser.parse_args()
input_dir = args.input_dir
sample_name = args.sample_name
threshold = float(args.threshold)
embeddings_file = args.embeddings_file
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
## 4. Load data
###############################################################################

# Load embeddings
with open(embeddings_file, 'rb') as f:
    dom_embeddings = pickle.load(f)
    
# Find files    
matching_files = []
for root, dirs, files in os.walk(input_dir):
    for file in files:
        if re.search(".*_annotdoms_resolved.tsv", file):
            matching_files.append(os.path.join(root, file))

# Load BGCs
bgc_domains = {}
for file in matching_files:
    bgc = os.path.basename(file).replace("_annotdoms_resolved.tsv", "")
    df_tmp = pd.read_csv(file, sep=' ', comment='#', header=None)
    bgc_domains[bgc] = df_tmp.iloc[:, 1].tolist()
    
###############################################################################
## 5. Map domains to embeddings 
###############################################################################

# Map domains to embeddings
bgc_embeddings = {}
for bgc in bgc_domains:
  
    # deduplicate elements in list
    bgc_domains[bgc] = list(dict.fromkeys(bgc_domains[bgc]))
    # initialize list of embeddings
    bgc_embeddings[bgc] = []
    
    for domain in bgc_domains[bgc]:
        if domain in dom_embeddings:
            bgc_embeddings[bgc].append(dom_embeddings[domain])
        # else:
        #     bgc_embeddings[bgc].append([0]*100) # if domain not in embeddings, add zeros

###############################################################################
## 6. Compute the mean embedding for each sample 
###############################################################################

bgc_embeddings_mean = {}
for bgc in bgc_embeddings:
    df_tmp = pd.DataFrame(bgc_embeddings[bgc]).T
    bgc_embeddings_mean[bgc] = df_tmp.mean(axis=1)
    
###############################################################################
## 7. Convert bgc_embeddings_mean to data frame
###############################################################################

df_bgc_embeddings_mean = pd.DataFrame(bgc_embeddings_mean).T

###############################################################################
# 8. Identify and remove rows with all NaNs
###############################################################################

i = df_bgc_embeddings_mean.isna().all(axis=1)
bgcs_ids_unannot = list(df_bgc_embeddings_mean[i].index)
df_bgc_embeddings_mean = df_bgc_embeddings_mean[~i]

###############################################################################
## 9. Cluster BGCs
###############################################################################

birch = Birch(
            n_clusters=None
        )
        
birch.threshold = threshold
birch.branching_factor = df_bgc_embeddings_mean.shape[0]

birch.fit(df_bgc_embeddings_mean)
brc_cluster = birch.predict(df_bgc_embeddings_mean)

###############################################################################
## 10. Format clustering as df
###############################################################################

bgc_ids = df_bgc_embeddings_mean.index

output_df = pd.DataFrame({
    'bgc_id': bgc_ids,
    'cluster_id': brc_cluster
})

###############################################################################
## 11. Add un annot BGCs
###############################################################################

cluster_id_max = output_df['cluster_id'].max()
cluster_ids_unannot = list(range(cluster_id_max +1, 
                                 cluster_id_max + len(bgcs_ids_unannot) +1))

output_unannot_df = pd.DataFrame({'bgc_id': bgcs_ids_unannot, 
                                  'cluster_id': str(cluster_ids_unannot)})

output_df = pd.concat([output_df, output_unannot_df], ignore_index=True)

###############################################################################
## 12. Export clustering
###############################################################################

output_df.to_csv(output_tsv, sep='\t', index=False)

print("cluster module executed successfully")
