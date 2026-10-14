#!/usr/bin/env python3
"""
SeDNA module 1: BGC annotation with antiSMASH.

- Runs antiSMASH (--minimal) on an assembly / MAG FASTA
- Moves the *.region*.gbk outputs into <output_dir>/gbks/
- Parses them into <output_dir>/annot_metadata.tsv (common schema shared by
  modules 1-3 and consumed by module 4)
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import (antismash_annot_parser, check_fasta, check_tools, find_files,
                   prepare_output_dir, run)

ANTISMASH_DB_DEFAULT = os.path.join(os.path.expanduser("~"), ".sedna", "db", "antismash")

###############################################################################
# 2. Define functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: BGC annotation with antiSMASH",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_fasta", dest="input_fasta", required=True,
        help="input assembly / MAG (FASTA)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name")

    parser.add_argument("--antismash_db", dest="antismash_db", default=ANTISMASH_DB_DEFAULT,
        help=f"antiSMASH databases directory (default: {ANTISMASH_DB_DEFAULT})")

    parser.add_argument("--taxon", dest="taxon", default="bacteria",
        choices=["bacteria", "fungi"],
        help="taxonomic classification of input sequences (default: bacteria)")

    parser.add_argument("--genefinding_tool", dest="genefinding_tool", default="prodigal-m",
        choices=["glimmerhmm", "prodigal", "prodigal-m", "none", "error"],
        help="antiSMASH gene finding algorithm (default: prodigal-m)")

    parser.add_argument("--minlength", dest="minlength", type=int, default=1000,
        help="only process sequences longer than this (default: 1000)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Run antiSMASH
###############################################################################

def run_antismash(input_fasta: str, output_dir: str, antismash_db: str, nslots: int,
                  genefinding_tool: str, taxon: str, minlength: int) -> None:
    run([
        "antismash",
        "--cpus", str(nslots),
        "--databases", antismash_db,
        "--genefinding-tool", genefinding_tool,
        "--taxon", taxon,
        "--allow-long-headers",
        "--minlength", str(minlength),
        "--minimal",
        "--output-dir", output_dir,
        input_fasta
    ])

###############################################################################
# 3. Main
###############################################################################

def main() -> None:

    args = parse_args()

    ###########################################################################
    # 3.1 Sanity checks
    ###########################################################################

    check_tools(["antismash"])
    check_fasta(args.input_fasta)

    if args.nslots <= 0:
        print("Error: --nslots must be a positive integer", file=sys.stderr)
        sys.exit(1)

    if args.minlength <= 0:
        print("Error: --minlength must be a positive integer", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.antismash_db):
        print(f"Error: antiSMASH database directory {args.antismash_db} not found. "
              "Run bin/download_dbs.py first", file=sys.stderr)
        sys.exit(1)

    ###########################################################################
    # 3.2 Create output dirs
    ###########################################################################

    prepare_output_dir(args.output_dir, args.overwrite)
    raw_output_dir = os.path.join(args.output_dir, "output")
    gbks_dir = os.path.join(args.output_dir, "gbks")
    os.makedirs(gbks_dir)

    ###########################################################################
    # 3.3 Run antiSMASH (it creates raw_output_dir itself and refuses an
    # existing one)
    ###########################################################################

    run_antismash(input_fasta=args.input_fasta,
                  output_dir=raw_output_dir,
                  antismash_db=args.antismash_db,
                  nslots=args.nslots,
                  genefinding_tool=args.genefinding_tool,
                  taxon=args.taxon,
                  minlength=args.minlength)

    ###########################################################################
    # 3.4 Move region GBKs to the gbks folder
    ###########################################################################

    antismash_output_gbks = find_files(input_dir=raw_output_dir, pattern=".*.region.*.gbk")

    if not antismash_output_gbks:
        print("No GBK files were generated by antiSMASH; no metadata file written")
        return

    for gbk in antismash_output_gbks:
        shutil.move(gbk, gbks_dir)

    ###########################################################################
    # 3.5 Parse and export metadata
    ###########################################################################

    antismash_annot_metadata = antismash_annot_parser(input_dir=gbks_dir,
                                                      sample_name=args.sample_name,
                                                      input_fasta=args.input_fasta)

    if antismash_annot_metadata is not None:
        antismash_annot_metadata.to_csv(os.path.join(args.output_dir, "annot_metadata.tsv"),
                                        sep='\t', index=False)
    else:
        print("No metadata file was generated after running antiSMASH")

    print(f"{os.path.basename(__file__)} executed successfully")


if __name__ == "__main__":
    main()
