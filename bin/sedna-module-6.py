#!/usr/bin/env python3
"""
SeDNA module 6: BGC embedding and clustering into Gene Cluster Families (GCFs).

Runs once over the module 5 outputs of every sample.

- Each BGC is the deduplicated list of Pfam domains in its
  *_annotdoms_resolved.tsv, mapped onto word2vec domain embeddings (MIBiG 3.1)
- --embedding_method sif (default, formerly cluster_sif): Smooth Inverse
  Frequency weighted sum, L2-normalised, with the first principal component
  removed (Arora et al., ICLR 2017)
- --embedding_method mean (formerly cluster_mean): unweighted mean, no
  normalisation, no PC removal
- BGCs are clustered with BIRCH (n_clusters=None: leaf subclusters are the
  GCFs; --threshold is the subcluster radius). BGCs with no mappable domain are
  appended as singleton GCFs
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import pickle
import re
import sys

import numpy as np
import pandas as pd
from sklearn.cluster import Birch
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(__file__))
from utils import check_file, find_files, prepare_output_dir

RESOURCES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                             "resources", "embeddings")
EMBEDDINGS_DEFAULT = os.path.join(RESOURCES_DIR, "dom_embeddings_mtx_mibig_pfam.pkl")
COUNTS_DEFAULT = os.path.join(RESOURCES_DIR, "domain2counts_mibig3.1_vs_pfam.pkl")

# SIF: domains occurring fewer times than this in MIBiG are dropped from the
# weight table
MIN_DOMAIN_COUNT = 10

###############################################################################
# 2. Define functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: BGC embedding and BIRCH clustering into GCFs",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_dirs", dest="input_dirs", nargs="+", required=True,
        help="module 5 output directories holding *_annotdoms_resolved.tsv files")

    parser.add_argument("--embedding_method", dest="embedding_method", default="sif",
        choices=["sif", "mean"],
        help="BGC embedding: sif (SIF-weighted, PC removed) or mean (unweighted) (default: sif)")

    parser.add_argument("--threshold", dest="threshold", type=float, default=1.0,
        help="BIRCH clustering threshold, i.e. GCF radius in embedding space (default: 1.0)")

    parser.add_argument("--alpha", dest="alpha", type=float, default=0.001,
        help="SIF a parameter (default: 0.001; --embedding_method sif only)")

    parser.add_argument("--embeddings_file", dest="embeddings_file", default=EMBEDDINGS_DEFAULT,
        help=f"domain embeddings pkl file (default: {EMBEDDINGS_DEFAULT})")

    parser.add_argument("--counts_file", dest="counts_file", default=COUNTS_DEFAULT,
        help=f"domain counts pkl file (default: {COUNTS_DEFAULT}; --embedding_method sif only)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Load the domain list of every BGC
###############################################################################

def load_bgc_domains(input_dirs: list) -> dict:
    """Return {bgc_id: [domain, ...]} from every non-empty *_annotdoms_resolved.tsv.

    BIRCH is sensitive to insertion order, so the BGCs are sorted by (sample
    name, file name). Otherwise the order would follow Nextflow's task
    completion order (the order of --input_dirs) and the GCFs could change
    from run to run on identical input.
    """
    files = []
    for input_dir in input_dirs:
        files.extend(find_files(input_dir=input_dir, pattern=".*_annotdoms_resolved.tsv"))
    # bgc_id = <tool>__<sample_name>__<contig>__<start>__<end>
    files.sort(key=lambda f: (os.path.basename(f).split("__")[1], os.path.basename(f)))

    bgc_domains = {}
    for file in files:
        if os.path.getsize(file) == 0:
            continue
        bgc = os.path.basename(file).replace("_annotdoms_resolved.tsv", "")
        df_tmp = pd.read_csv(file, sep=' ', comment='#', header=None)
        bgc_domains[bgc] = df_tmp.iloc[:, 1].tolist()
    return bgc_domains

###############################################################################
# 2.3 SIF: weights from domain counts
###############################################################################

def compute_sif_weights(counts_file: str, alpha: float) -> dict:
    with open(counts_file, 'rb') as f:
        dom2counts = pickle.load(f)

    dom2counts_dict = dom2counts.set_index('domain_id')['count'].to_dict()
    dom2counts_dict = {dom: count for dom, count in dom2counts_dict.items()
                       if count >= MIN_DOMAIN_COUNT}

    total_words = sum(dom2counts_dict.values())
    return {word: alpha / (alpha + (freq / total_words)) for word, freq in dom2counts_dict.items()}

###############################################################################
# 2.4 SIF: weighted BGC embedding
###############################################################################

def map_domain2embeddings(domains: list, domain_embeddings: dict, weights: dict) -> np.ndarray:
    """SIF-weighted sum of the domain vectors, L2-normalised (all-NaN if empty)."""
    vector_size = len(next(iter(domain_embeddings.values())))
    embedding_mean = np.zeros(vector_size)

    for domain in domains:
        if domain in domain_embeddings:
            weight = weights.get(domain, 1e-3)
            embedding_mean += weight * np.array(domain_embeddings[domain])

    embedding_norm = np.linalg.norm(embedding_mean)
    if embedding_norm == 0:
        return np.full(vector_size, np.nan)
    return embedding_mean / embedding_norm

###############################################################################
# 2.5 SIF: remove the first principal component
###############################################################################

def remove_pc(embeddings: np.ndarray) -> np.ndarray:
    pca = PCA(n_components=1)
    pca.fit(embeddings)
    pc = pca.components_
    return embeddings - embeddings.dot(pc.T) * pc

###############################################################################
# 2.6 SIF embeddings (NaN rows split off, PC removed)
###############################################################################

def embed_sif(bgc_domains: dict, domain_embeddings: dict, weights: dict):
    """Return (embeddings DF, unannotated BGC ids)."""
    bgc_embeddings = {}
    for bgc, domains in bgc_domains.items():
        # Deduplicate elements in list
        bgc_embeddings[bgc] = map_domain2embeddings(list(dict.fromkeys(domains)),
                                                    domain_embeddings, weights)
    embedding_df = pd.DataFrame(bgc_embeddings).T

    # Identify and remove rows with all NaNs
    i = embedding_df.isna().all(axis=1)
    bgcs_ids_unannot = list(embedding_df[i].index)
    embedding_df = embedding_df[~i]

    embedding_df = pd.DataFrame(remove_pc(embedding_df.to_numpy()),
                                index=embedding_df.index,
                                columns=embedding_df.columns)
    return embedding_df, bgcs_ids_unannot

###############################################################################
# 2.7 Mean embeddings (exported before NaN filtering, as cluster_mean did)
###############################################################################

def embed_mean(bgc_domains: dict, domain_embeddings: dict):
    """Return (embeddings DF incl. NaN rows, unannotated BGC ids)."""
    bgc_embeddings_mean = {}
    for bgc, domains in bgc_domains.items():
        vectors = [domain_embeddings[d] for d in dict.fromkeys(domains) if d in domain_embeddings]
        bgc_embeddings_mean[bgc] = pd.DataFrame(vectors).T.mean(axis=1)
    embedding_df = pd.DataFrame(bgc_embeddings_mean).T

    i = embedding_df.isna().all(axis=1)
    return embedding_df, list(embedding_df[i].index)

###############################################################################
# 2.8 BIRCH clustering, unannotated BGCs appended as singletons
###############################################################################

def cluster_birch(embedding_df: pd.DataFrame, bgcs_ids_unannot: list,
                  threshold: float, embedding_method: str) -> pd.DataFrame:
    if embedding_method == "sif":
        birch = Birch(n_clusters=None, compute_labels=False, copy=False)
    else:
        birch = Birch(n_clusters=None)

    birch.threshold = threshold
    birch.branching_factor = embedding_df.shape[0]
    birch.fit(embedding_df)
    brc_cluster = birch.predict(embedding_df)

    output_df = pd.DataFrame({'bgc_id': embedding_df.index, 'cluster_id': brc_cluster})

    cluster_id_max = output_df['cluster_id'].max()
    cluster_ids_unannot = [str(c) for c in range(cluster_id_max + 1,
                                                 cluster_id_max + len(bgcs_ids_unannot) + 1)]
    output_unannot_df = pd.DataFrame({'bgc_id': bgcs_ids_unannot,
                                      'cluster_id': cluster_ids_unannot})

    return pd.concat([output_df, output_unannot_df], ignore_index=True)

###############################################################################
# 3. Main
###############################################################################

def main() -> None:

    args = parse_args()

    ###########################################################################
    # 3.1 Sanity checks
    ###########################################################################

    check_file(args.embeddings_file, f"--embeddings_file {args.embeddings_file}")
    if args.embedding_method == "sif":
        check_file(args.counts_file, f"--counts_file {args.counts_file}")

    ###########################################################################
    # 3.2 Create output dir
    ###########################################################################

    prepare_output_dir(args.output_dir, args.overwrite)
    thres_str = re.sub(r'\.?0+$', '', f"{args.threshold:f}")
    bgc_embeddings_tsv = os.path.join(args.output_dir,
                                      f"bgc_embeddings-{args.embedding_method}.tsv")
    output_tsv = os.path.join(args.output_dir,
                              f"bgc_clust-{args.embedding_method}-thres{thres_str}.tsv")

    ###########################################################################
    # 3.3 Load data
    ###########################################################################

    with open(args.embeddings_file, 'rb') as f:
        domain_embeddings = pickle.load(f)

    bgc_domains = load_bgc_domains(args.input_dirs)
    if len(bgc_domains) == 0:
        print("Error: no non-empty *_annotdoms_resolved.tsv file found in --input_dirs",
              file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.4 Compute BGC embeddings
    ###########################################################################

    if args.embedding_method == "sif":
        weights = compute_sif_weights(args.counts_file, args.alpha)
        embedding_df, bgcs_ids_unannot = embed_sif(bgc_domains, domain_embeddings, weights)
        embedding_df.to_csv(bgc_embeddings_tsv, sep='\t', index=True)
    else:
        embedding_df, bgcs_ids_unannot = embed_mean(bgc_domains, domain_embeddings)
        embedding_df.to_csv(bgc_embeddings_tsv, sep='\t', index=True)
        embedding_df = embedding_df.drop(index=bgcs_ids_unannot)

    ###########################################################################
    # 3.5 Cluster BGCs and export
    ###########################################################################

    output_df = cluster_birch(embedding_df, bgcs_ids_unannot,
                              args.threshold, args.embedding_method)
    output_df.to_csv(output_tsv, sep='\t', index=False)

    print(f"{os.path.basename(__file__)} executed successfully")


if __name__ == "__main__":
    main()
