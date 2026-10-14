#!/usr/bin/env python3
"""
SeDNA module 4: Dereplication of BGCs predicted by different tools.

- Reads the annot_metadata.tsv of each module 1-3 output directory given in
  --input_dirs (missing tables, i.e. tools that found nothing, are skipped)
- Merges the tables pairwise and recursively with coordinate-based
  dereplication (13 topological cases, longest-wins, --overlap_thres)
- Writes four sibling directories, each with gbks/ and annot_metadata.tsv:
  dereplicated/ (the catalog), overlapped/, partially_overlapped/,
  non_overlapped/ (audit)
- Adds bgc_class_formatted to the dereplicated catalog, mapping each tool's
  class vocabulary through --synonyms_file
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from utils import (check_file, copy_gbks, create_df, dereplicate_non_shared_contigs,
                   dereplicate_shared_contigs, map_synonyms, prepare_output_dir,
                   sort_single_table)

SYNONYMS_DEFAULT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "resources", "bgc_class_synonyms.tsv")

OUTPUT_SUBDIRS = ["dereplicated", "overlapped", "partially_overlapped", "non_overlapped"]

###############################################################################
# 2. Define functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: Dereplicate BGCs annotated with different tools",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_dirs", dest="input_dirs", nargs="+", required=True,
        help="module 1-3 output directories, each holding gbks/ and (if any BGC "
             "was found) annot_metadata.tsv. Merged in the order given")

    parser.add_argument("--overlap_thres", dest="overlap_thres", type=float, default=0.75,
        help="fraction of overlap (relative to the longest BGC) above which two "
             "staggered BGCs are called overlapped rather than partially overlapped "
             "(default: 0.75)")

    parser.add_argument("--synonyms_file", dest="synonyms_file", default=SYNONYMS_DEFAULT,
        help=f"table of standardized BGC class names across tools (default: {SYNONYMS_DEFAULT})")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Load a tool's metadata table, pointing `file` at the staged GBKs
###############################################################################

def load_metadata(input_dir: str) -> pd.DataFrame:
    """Return the tool's metadata DF, or None if the tool found no BGC.

    The `file` column holds the absolute path the GBK had inside the upstream
    task's work dir, which does not exist in this task. The GBK is staged here
    under <input_dir>/gbks/ with the same basename, so the column is remapped.
    """
    metadata_tsv = os.path.join(input_dir, "annot_metadata.tsv")
    if not os.path.isfile(metadata_tsv):
        print(f"No annot_metadata.tsv in {input_dir}; skipping")
        return None

    metadata_df = pd.read_csv(metadata_tsv, sep='\t')
    if metadata_df.shape[0] == 0:
        return None

    gbks_dir = os.path.join(os.path.abspath(input_dir), "gbks")
    metadata_df['file'] = [os.path.join(gbks_dir, os.path.basename(f)) for f in metadata_df['file']]
    return metadata_df

###############################################################################
# 2.3 Recursive pairwise dereplication
###############################################################################

def recursive_dereplication(metadata: list,
                            dereplicated_bgcs: dict,
                            overlapped_bgcs: dict,
                            partially_overlapped_bgcs: dict,
                            non_overlapped_bgcs: dict,
                            overlap_thres: float,
                            output_dir: str) -> dict:
    """Merge metadata[0] with metadata[1]; the merged catalog becomes the new
    metadata[0], and recurse until a single table is left."""

    metadata_df1 = metadata[0]
    metadata_df2 = metadata[1]

    # dereplicate shared contigs
    output_shared_contigs = dereplicate_shared_contigs(metadata_df1=metadata_df1,
                                                       metadata_df2=metadata_df2,
                                                       dereplicated_bgcs=dereplicated_bgcs,
                                                       overlapped_bgcs=overlapped_bgcs,
                                                       partially_overlapped_bgcs=partially_overlapped_bgcs,
                                                       non_overlapped_bgcs=non_overlapped_bgcs,
                                                       overlap_thres=overlap_thres)

    dereplicated_bgcs = output_shared_contigs['dereplicated_bgcs']
    overlapped_bgcs = output_shared_contigs['overlapped_bgcs']
    partially_overlapped_bgcs = output_shared_contigs['partially_overlapped_bgcs']
    non_overlapped_bgcs = output_shared_contigs['non_overlapped_bgcs']

    # dereplicate non-shared contigs
    output_non_shared_contigs = dereplicate_non_shared_contigs(metadata_df1=metadata_df1,
                                                               metadata_df2=metadata_df2,
                                                               dereplicated_bgcs=dereplicated_bgcs,
                                                               non_overlapped_bgcs=non_overlapped_bgcs)

    dereplicated_bgcs = output_non_shared_contigs['dereplicated_bgcs']
    non_overlapped_bgcs = output_non_shared_contigs['non_overlapped_bgcs']

    # update metadata list: the merged catalog replaces the first two tables
    metadata = metadata[1:]
    metadata[0] = create_df(input_dict=dereplicated_bgcs,
                            output_dir=os.path.join(output_dir, "dereplicated"),
                            output_dir_gbks=os.path.join(output_dir, "dereplicated", "gbks"),
                            df=True)

    if len(metadata) == 1:
        return {'dereplicated_bgcs': dereplicated_bgcs,
                'overlapped_bgcs': overlapped_bgcs,
                'partially_overlapped_bgcs': partially_overlapped_bgcs,
                'non_overlapped_bgcs': non_overlapped_bgcs}

    return recursive_dereplication(metadata=metadata,
                                   dereplicated_bgcs=dereplicated_bgcs,
                                   overlapped_bgcs=overlapped_bgcs,
                                   partially_overlapped_bgcs=partially_overlapped_bgcs,
                                   non_overlapped_bgcs=non_overlapped_bgcs,
                                   overlap_thres=overlap_thres,
                                   output_dir=output_dir)

###############################################################################
# 3. Main
###############################################################################

def main() -> None:

    args = parse_args()

    ###########################################################################
    # 3.1 Sanity checks
    ###########################################################################

    check_file(args.synonyms_file, f"--synonyms_file {args.synonyms_file}")
    for input_dir in args.input_dirs:
        if not os.path.isdir(input_dir):
            print(f"Error: input dir {input_dir} does not exist", file=sys.stderr)
            sys.exit(1)

    ###########################################################################
    # 3.2 Create output dirs
    ###########################################################################

    prepare_output_dir(args.output_dir, args.overwrite)
    for subdir in OUTPUT_SUBDIRS:
        os.makedirs(os.path.join(args.output_dir, subdir, "gbks"))

    ###########################################################################
    # 3.3 Load the metadata tables
    ###########################################################################

    metadata_list = []
    for input_dir in args.input_dirs:
        metadata_df = load_metadata(input_dir)
        if metadata_df is not None:
            metadata_list.append(metadata_df)

    if len(metadata_list) == 0:
        # Not an error: a sample can legitimately carry no BGC. The empty
        # catalog keeps the directory layout downstream modules expect.
        print("No BGCs found by any tool; writing an empty catalog")
        return

    ###########################################################################
    # 3.4 Initialize dictionaries
    ###########################################################################

    tools = [metadata_df['tool'].iloc[0] for metadata_df in metadata_list]
    dereplicated_bgcs = {tool: dict() for tool in tools}
    overlapped_bgcs = {tool: dict() for tool in tools}
    partially_overlapped_bgcs = {tool: dict() for tool in tools}
    non_overlapped_bgcs = {tool: dict() for tool in tools}

    ###########################################################################
    # 3.5 Execute dereplication
    ###########################################################################

    if len(metadata_list) > 1:
        output = recursive_dereplication(metadata=metadata_list,
                                         dereplicated_bgcs=dereplicated_bgcs,
                                         overlapped_bgcs=overlapped_bgcs,
                                         partially_overlapped_bgcs=partially_overlapped_bgcs,
                                         non_overlapped_bgcs=non_overlapped_bgcs,
                                         overlap_thres=args.overlap_thres,
                                         output_dir=args.output_dir)
    else:
        sort_single_table(metadata_df1=metadata_list[0],
                          dereplicated_bgcs=dereplicated_bgcs,
                          non_overlapped_bgcs=non_overlapped_bgcs)
        output = {'dereplicated_bgcs': dereplicated_bgcs,
                  'overlapped_bgcs': overlapped_bgcs,
                  'partially_overlapped_bgcs': partially_overlapped_bgcs,
                  'non_overlapped_bgcs': non_overlapped_bgcs}

    ###########################################################################
    # 3.6 Sort data: copy GBKs
    ###########################################################################

    output_keys = {"dereplicated": "dereplicated_bgcs",
                   "overlapped": "overlapped_bgcs",
                   "partially_overlapped": "partially_overlapped_bgcs",
                   "non_overlapped": "non_overlapped_bgcs"}

    for subdir, key in output_keys.items():
        copy_gbks(input_dict=output[key],
                  output_dir=os.path.join(args.output_dir, subdir, "gbks"))

    ###########################################################################
    # 3.7 Map formatted BGC class names
    ###########################################################################

    output['dereplicated_bgcs'] = map_synonyms(synonyms_file=args.synonyms_file,
                                               input_dict=output['dereplicated_bgcs'])

    ###########################################################################
    # 3.8 Write metadata tables
    ###########################################################################

    for subdir, key in output_keys.items():
        create_df(input_dict=output[key],
                  output_dir=os.path.join(args.output_dir, subdir),
                  output_dir_gbks=os.path.join(args.output_dir, subdir, "gbks"))

    print(f"{os.path.basename(__file__)} executed successfully")


if __name__ == "__main__":
    main()
