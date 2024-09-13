#!/usr/bin/env python3

###############################################################################
## 1. Set env
###############################################################################

import argparse
import pandas
import re
import os
import pickle
from sklearn.cluster import Birch

###############################################################################
## 2. Parse input data and optional arguments
###############################################################################

parser = argparse.ArgumentParser(prog='cluster.py', \
                                 description='Cluster BGCs into GCFs.')

# general parameters
parser.add_argument("--input_sample", help="Input fasta file.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--output_dir", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")

# Get general parameters
args = parser.parse_args()
input_sample = args.input_sample
output_dir = args.output_dir
sample_name = args.sample_name

###############################################################################
## 3. Load data
###############################################################################

# Load embeddings
with open('results/embeddings/dom_embeddings_mtx_mibig_pfam.pkl', 'rb') as f:
    dom_embeddings = pickle.load(f)
    
# find annotdoms_resolved.tsv files
input_dir = "/home/epereira/workspace/dev/new_atlantis/dev/clustering/data/export_testdata/"

matching_files = []
for root, dirs, files in os.walk(input_dir):
    for file in files:
        if re.search(".*_annotdoms_resolved.tsv", file):
            matching_files.append(os.path.join(root, file))


# Load BGCs
bgc_domains = {}
for file in matching_files:
    bgc = os.path.basename(file).replace("_annotdoms_resolved.tsv", "")
    df_tmp = pandas.read_csv(file, sep=' ', comment='#', header=None)
    bgc_domains[bgc] = df_tmp.iloc[:, 1].tolist()

###############################################################################
## 4. Map domains to embeddings 
###############################################################################

# Map domains to embeddings
bgc_embeddings = {}
for bgc in bgc_domains:
  
    # deduplicate ements in list
    bgc_domains[bgc] = list(dict.fromkeys(bgc_domains[bgc]))
    # initialize list of embeddings
    bgc_embeddings[bgc] = []
    
    for domain in bgc_domains[bgc]:
        if domain in dom_embeddings:
            bgc_embeddings[bgc].append(dom_embeddings[domain])
        # else:
        #     bgc_embeddings[bgc].append([0]*100) # if domain not in embeddings, add zeros

###############################################################################
## 5. Compute the mean embedding for each sample 
###############################################################################

bgc_embeddings_mean = {}
for bgc in bgc_embeddings:
    df_tmp = pd.DataFrame(bgc_embeddings[bgc]).T
    bgc_embeddings_mean[bgc] = df_tmp.mean(axis=1)
    
###############################################################################
## 6. Convert bgc_embeddings_mean to data frame
###############################################################################

df_bgc_embeddings_mean = pd.DataFrame(bgc_embeddings_mean).T

###############################################################################
## 7. Cluster BGCs
###############################################################################

brc = Birch(branching_factor=5, n_clusters=None, 
            threshold=0.5, compute_labels=True)
brc.fit(df_bgc_embeddings_mean)
