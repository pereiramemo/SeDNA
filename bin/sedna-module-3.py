#!/usr/bin/env python3
"""
SeDNA module 3: BGC annotation with GECCO.

- Runs `gecco run` on an assembly / MAG FASTA
- Copies the *_cluster_<n>.gbk outputs into <output_dir>/gbks/, renamed to the
  antiSMASH convention (<contig>.region<NNN>.gbk)
- Parses GECCO's *.clusters.tsv into <output_dir>/annot_metadata.tsv
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import (check_fasta, check_tools, find_files, gecco_annot_parser,
                   prepare_output_dir, run)

###############################################################################
# 2. Define functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: BGC annotation with GECCO",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_fasta", dest="input_fasta", required=True,
        help="input assembly / MAG (FASTA)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name")

    parser.add_argument("--score_thres", dest="score_thres", type=float, default=0.80,
        help="probability threshold for cluster detection (default: 0.80)")

    parser.add_argument("--cds_count_thres", dest="cds_count_thres", type=int, default=3,
        help="minimum number of CDSs a valid cluster must contain (default: 3)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Run GECCO
###############################################################################

def run_gecco(input_fasta: str, output_dir: str, nslots: int,
              cds_count_thres: int, score_thres: float) -> None:
    run([
        "gecco", "run",
        "--genome", input_fasta,
        "--cds", str(cds_count_thres),
        "--threshold", str(score_thres),
        "--jobs", str(nslots),
        "--output-dir", output_dir
    ])

###############################################################################
# 2.3 Copy GECCO GBKs, renamed to the antiSMASH region convention
###############################################################################

def move_gbks(gecco_gbk_list: list, output_dir: str) -> None:
    for file_path in gecco_gbk_list:
        file = os.path.basename(file_path)
        match = re.search(r'_cluster_(\d+)\.gbk', file)
        number = int(match.group(1))
        file_renamed = re.sub(r'_cluster_\d+.gbk', '.region{:03d}.gbk'.format(number), file)
        shutil.copy(file_path, os.path.join(output_dir, file_renamed))

###############################################################################
# 3. Main
###############################################################################

def main() -> None:

    args = parse_args()

    ###########################################################################
    # 3.1 Sanity checks
    ###########################################################################

    check_tools(["gecco"])
    check_fasta(args.input_fasta)

    if args.nslots <= 0:
        print("Error: --nslots must be a positive integer", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.2 Create output dirs
    ###########################################################################

    prepare_output_dir(args.output_dir, args.overwrite)
    raw_output_dir = os.path.join(args.output_dir, "output")
    gbks_dir = os.path.join(args.output_dir, "gbks")
    os.makedirs(raw_output_dir)
    os.makedirs(gbks_dir)

    ###########################################################################
    # 3.3 Run GECCO
    ###########################################################################

    run_gecco(input_fasta=args.input_fasta,
              output_dir=raw_output_dir,
              nslots=args.nslots,
              cds_count_thres=args.cds_count_thres,
              score_thres=args.score_thres)

    ###########################################################################
    # 3.4 Copy GBKs to the gbks folder, renamed as *.region*.gbk
    ###########################################################################

    gecco_gbk_list = find_files(input_dir=raw_output_dir, pattern=r'.cluster_[0-9]+.gbk')
    move_gbks(gecco_gbk_list=gecco_gbk_list, output_dir=gbks_dir)

    ###########################################################################
    # 3.5 Parse and export metadata. GECCO's own clusters.tsv is the source,
    # so the parser reads raw_output_dir, not the gbks folder.
    ###########################################################################

    gecco_annot_metadata = gecco_annot_parser(input_dir=raw_output_dir,
                                              sample_name=args.sample_name,
                                              input_fasta=args.input_fasta)

    if gecco_annot_metadata is not None:
        gecco_annot_metadata.to_csv(os.path.join(args.output_dir, "annot_metadata.tsv"),
                                    sep='\t', index=False)
    else:
        print("No metadata file was generated after running GECCO")

    print(f"{os.path.basename(__file__)} executed successfully")


if __name__ == "__main__":
    main()
