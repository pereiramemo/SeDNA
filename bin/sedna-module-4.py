#!/usr/bin/env python3
"""
SeDNA module 4: Dereplication of BGCs predicted by different tools.

- Reads the annot_metadata.tsv of each module 1-3 output directory given in
  --input_dirs (missing tables, i.e. tools that found nothing, are skipped)
- Merges the tables pairwise and recursively with coordinate-based
  dereplication (13 topological cases, longest-wins, --overlap_thres)
- Writes the catalog in three versions, each with gbks/ and annot_metadata.tsv:
  dereplicated/longest/  the longest prediction of each locus (the result of
                         the pairwise dereplication)
  dereplicated/core/     the region supported by the most tools (the
                         intersection when all tools overlap), longest on ties
  dereplicated/extended/ leftmost start to rightmost end of all supporting
                         predictions
  plus the audit dirs overlapped/, partially_overlapped/, non_overlapped/
- Adds bgc_class_formatted to the catalog, mapping each tool's class
  vocabulary through --synonyms_file
- Works for any number of tools: core and extended are derived from the
  longest version and the input tables only, never from tool names
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import gzip
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from utils import (check_file, collect_cds, copy_gbks, create_df, dereplicate_non_shared_contigs,
                   dereplicate_shared_contigs, find_bgc_groups, is_gzipped, make_bgc_id,
                   map_synonyms, most_supported_region, prepare_output_dir, sort_single_table,
                   write_region_gbk)

SYNONYMS_DEFAULT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "resources", "bgc_class_synonyms.tsv")

DEREP_VERSIONS = ["longest", "core", "extended"]
OUTPUT_SUBDIRS = [os.path.join("dereplicated", v) for v in DEREP_VERSIONS] + \
                 ["overlapped", "partially_overlapped", "non_overlapped"]

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

    parser.add_argument("--input_fasta", dest="input_fasta", required=True,
        help="assembly FASTA of the sample (optionally gzipped); the core and "
             "extended GBKs are cut from it")

    parser.add_argument("--overlap_thres", dest="overlap_thres", type=float, default=0.80,
        help="fraction of overlap (relative to the longest BGC) above which two "
             "staggered BGCs are called overlapped rather than partially overlapped "
             "(default: 0.80)")

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

    # Ids are read as text: otherwise contigs named 1, 2 or NA become numbers
    # or NaN and are not found among the FASTA ids
    metadata_df = pd.read_csv(metadata_tsv, sep='\t',
                              dtype={'acc': str, 'contig_id': str, 'sample_name': str},
                              keep_default_na=False)
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
                            output_dir=os.path.join(output_dir, "dereplicated", "longest"),
                            output_dir_gbks=os.path.join(output_dir, "dereplicated", "longest", "gbks"),
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
# 2.4 Load the contig sequences of the sample
###############################################################################

def load_contigs(input_fasta: str) -> dict:
    from Bio import SeqIO
    opener = gzip.open if is_gzipped(input_fasta) else open
    with opener(input_fasta, "rt") as handle:
        return {record.id: str(record.seq) for record in SeqIO.parse(handle, "fasta")}

###############################################################################
# 2.5 Build the longest, core and extended versions of the catalog
###############################################################################

def build_versions(dereplicated_bgcs: dict, metadata_list: list, contigs: dict,
                   output_dir: str) -> dict:
    """Return {version: {key: {idx: {'row': Series}}}} for create_df, after
    writing the core and extended GBKs.

    The groups come from the longest version and all input predictions only,
    so this works for any number of tools. Every version gets `tools` (the
    supporting tools, in input order, joined by '+'), `n_tools` and `members`;
    core also gets `core_support`. core and extended take tool, bgc_class and
    bgc_class_formatted from the representative and use `tools` as bgc_id
    prefix, so a BGC predicted by one tool has the same bgc_id in all three.
    """

    # All input predictions, in input-table order
    preds = []
    tool_rank = {}
    for metadata_df in metadata_list:
        for _, row in metadata_df.iterrows():
            row_as_dict = row.to_dict()
            tool_rank.setdefault(row_as_dict['tool'], len(tool_rank))
            preds.append({'bgc_id': make_bgc_id(row_as_dict),
                          'contig_id': row_as_dict['contig_id'],
                          'start': int(row_as_dict['start']),
                          'end': int(row_as_dict['end']),
                          'tool': row_as_dict['tool'],
                          'file': row_as_dict['file']})
    preds_by_id = {p['bgc_id']: p for p in preds}

    # Representatives (the longest version)
    reps = []
    for key in dereplicated_bgcs:
        for idx in dereplicated_bgcs[key]:
            row = dereplicated_bgcs[key][idx]['row'].copy()
            reps.append((key, idx, make_bgc_id(row.to_dict()), row))

    groups = find_bgc_groups(preds, [rep_id for _, _, rep_id, _ in reps])

    versions = {v: {} for v in ["longest", "core", "extended"]}
    for key, idx, rep_id, rep_row in reps:
        members = [preds_by_id[m] for m in groups[rep_id]]
        tools = sorted(set(m['tool'] for m in members), key=lambda t: tool_rank[t])

        rep_row['tools'] = '+'.join(tools)
        rep_row['n_tools'] = len(tools)
        rep_row['members'] = ','.join(groups[rep_id])
        versions["longest"].setdefault(key, {})[idx] = {'row': rep_row}

        contig_id = rep_row['contig_id']
        contig_length = int(rep_row['contig_length'])
        contig_seq = contigs[contig_id]

        cds_features = None
        if len(members) > 1:
            cds_features = collect_cds([(m['file'], m['start']) for m in members], contig_seq)

        core_start, core_end, core_support = most_supported_region(
            [(m['start'], m['end'], m['tool']) for m in members])
        coords = {"core": (core_start, core_end),
                  "extended": (min(m['start'] for m in members), max(m['end'] for m in members))}

        for version, (start, end) in coords.items():
            row = rep_row.copy()
            row['start'] = start
            row['end'] = end
            row['length'] = end - start
            if len(members) > 1:
                row['on_edge'] = str(start < 3 or (contig_length - end) <= 30)
            if version == "core":
                row['core_support'] = core_support

            bgc_id = make_bgc_id(row.to_dict(), 'tools')
            gbk = os.path.join(output_dir, "dereplicated", version, "gbks", f"{bgc_id}.gbk")
            if os.path.exists(gbk):
                # two representatives sharing members can yield the same region
                print(f"Warning: {version} BGC {bgc_id} produced twice; keeping one", file=sys.stderr)
                continue
            if len(members) > 1:
                write_region_gbk(contig_id=contig_id, contig_seq=contig_seq, start=start, end=end,
                                 cds_features=cds_features, bgc_class=row['bgc_class'],
                                 on_edge=row['on_edge'], tool=row['tools'], output_gbk=gbk)
            else:
                # one prediction: the region is the prediction itself
                shutil.copyfile(members[0]['file'], gbk)
            versions[version].setdefault(key, {})[idx] = {'row': row}

    return versions

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

    output_keys = {os.path.join("dereplicated", "longest"): "dereplicated_bgcs",
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
    # 3.8 Build the core and extended versions
    ###########################################################################

    versions = build_versions(dereplicated_bgcs=output['dereplicated_bgcs'],
                              metadata_list=metadata_list,
                              contigs=load_contigs(args.input_fasta),
                              output_dir=args.output_dir)

    ###########################################################################
    # 3.9 Write metadata tables
    ###########################################################################

    for subdir, key in output_keys.items():
        if key == "dereplicated_bgcs":
            continue
        create_df(input_dict=output[key],
                  output_dir=os.path.join(args.output_dir, subdir),
                  output_dir_gbks=os.path.join(args.output_dir, subdir, "gbks"))

    for version in DEREP_VERSIONS:
        version_dir = os.path.join(args.output_dir, "dereplicated", version)
        create_df(input_dict=versions[version],
                  output_dir=version_dir,
                  output_dir_gbks=os.path.join(version_dir, "gbks"),
                  id_prefix_field="tool" if version == "longest" else "tools")

    print(f"{os.path.basename(__file__)} executed successfully")


if __name__ == "__main__":
    main()
