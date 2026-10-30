#!/usr/bin/env python3
"""
SeDNA module 2: BGC annotation with deepBGC.

Runs inside the module-2 image, which is Python 3.7 (deepBGC needs
TensorFlow 1.15); keep this script and utils.py free of newer syntax.

- Runs `deepbgc pipeline` (deepbgc detector + product_class classifier)
- Splits its single multi-record GBK into one <contig>.region<NNN>.gbk per BGC
- Filters low-confidence BGCs post hoc (--score_thres, --cds_count_thres) into
  gbks_removed/, so stringency can be relaxed later without re-running deepBGC
- Parses the kept GBKs into <output_dir>/annot_metadata.tsv and rewrites them
  as antiSMASH-like region GBKs
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import (check_fasta, check_tools, deepbgc_annot_parser, find_files,
                   format_gbks, gbk_filter, gbk_splitter, prepare_output_dir, run)

DEEPBGC_DB_DEFAULT = os.path.join(os.path.expanduser("~"), ".sedna", "db", "deepbgc", "data")

###############################################################################
# 2. Define functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: BGC annotation with deepBGC",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_fasta", dest="input_fasta", required=True,
        help="input assembly / MAG (FASTA)")

    parser.add_argument("--sample_name", dest="sample_name", required=True,
        help="sample name")

    parser.add_argument("--deepbgc_db", dest="deepbgc_db", default=DEEPBGC_DB_DEFAULT,
        help=f"deepBGC downloads directory, exported as DEEPBGC_DOWNLOADS_DIR "
             f"(default: {DEEPBGC_DB_DEFAULT})")

    parser.add_argument("--score_thres", dest="score_thres", type=float, default=0.90,
        help="minimum deepBGC score to keep a BGC (default: 0.90)")

    parser.add_argument("--cds_count_thres", dest="cds_count_thres", type=int, default=3,
        help="minimum number of CDSs to keep a BGC (default: 3)")

    parser.add_argument("--minlength", dest="minlength", type=int, default=1000,
        help="minimum BGC nucleotide length (default: 1000)")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads (default: 4; deepBGC itself is single-threaded)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 Run deepBGC
###############################################################################

def run_deepbgc(input_fasta: str, output_dir: str, minlength: int) -> None:
    run([
        "deepbgc", "pipeline",
        "--detector", "deepbgc",
        "--classifier", "product_class",
        "--prodigal-meta-mode",
        "--min-nucl", str(minlength),
        "--output", output_dir,
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

    check_tools(["deepbgc"])
    check_fasta(args.input_fasta)

    if args.minlength <= 0:
        print("Error: --minlength must be a positive integer", file=sys.stderr)
        sys.exit(1)

    if not os.path.isdir(args.deepbgc_db):
        print(f"Error: deepBGC downloads directory {args.deepbgc_db} not found. "
              "Run bin/download_dbs.py first", file=sys.stderr)
        sys.exit(1)

    # deepBGC looks up its models and Pfam DB here instead of ~/.local/share/deepbgc
    os.environ["DEEPBGC_DOWNLOADS_DIR"] = args.deepbgc_db

    ###########################################################################
    # 3.2 Create output dirs
    ###########################################################################

    prepare_output_dir(args.output_dir, args.overwrite)
    raw_output_dir = os.path.join(args.output_dir, "output")
    gbks_dir = os.path.join(args.output_dir, "gbks")
    gbks_removed_dir = os.path.join(args.output_dir, "gbks_removed")
    os.makedirs(raw_output_dir)

    ###########################################################################
    # 3.3 Run deepBGC
    ###########################################################################

    run_deepbgc(input_fasta=args.input_fasta,
                output_dir=raw_output_dir,
                minlength=args.minlength)

    ###########################################################################
    # 3.4 Find the (single) GBK output file and split it per BGC
    ###########################################################################

    deepbgc_output_gbk = find_files(input_dir=raw_output_dir, pattern=".bgc.gbk")

    if len(deepbgc_output_gbk) > 1:
        print("Error: more than one deepBGC GBK output file found; there should be only one:\n"
              + "\n".join(deepbgc_output_gbk), file=sys.stderr)
        sys.exit(1)

    if len(deepbgc_output_gbk) == 0:
        print("No GBK files were generated by deepBGC; no metadata file written")
        return

    gbk_splitter(input_file=deepbgc_output_gbk[0],
                 output_dir=gbks_dir,
                 sample_name=args.sample_name)

    ###########################################################################
    # 3.5 Filter deepBGC outputs with low score or few CDSs
    ###########################################################################

    gbk_filter(input_dir=gbks_dir,
               sample_name=args.sample_name,
               deepbgc_score_thres=args.score_thres,
               deepbgc_cds_count_thres=args.cds_count_thres,
               output_dir=gbks_removed_dir)

    ###########################################################################
    # 3.6 Parse and export metadata
    ###########################################################################

    deepbgc_annot_metadata = deepbgc_annot_parser(input_dir=gbks_dir,
                                                  sample_name=args.sample_name,
                                                  input_fasta=args.input_fasta)

    if deepbgc_annot_metadata is None:
        print("No metadata file was generated after running deepBGC")
        return

    deepbgc_annot_metadata_tsv = os.path.join(args.output_dir, "annot_metadata.tsv")
    deepbgc_annot_metadata.to_csv(deepbgc_annot_metadata_tsv, sep='\t', index=False)

    ###########################################################################
    # 3.7 Format deepBGC GBKs as antiSMASH-like region GBKs
    ###########################################################################

    format_gbks(input_dir=gbks_dir,
                input_tsv=deepbgc_annot_metadata_tsv,
                output_dir=gbks_dir,
                tool="deepbgc")

    print(f"{os.path.basename(__file__)} executed successfully")


if __name__ == "__main__":
    main()
