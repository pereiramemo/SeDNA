#!/usr/bin/env python3
"""
SeDNA module 5: Protein-domain annotation of the dereplicated BGC catalog.

Two modes, writing the same <bgc_id>_annotdoms_resolved.tsv files (one per BGC)
that module 6 consumes:

- default (annot_cds): extract the translated CDSs of every GBK in --input_dir,
  search them against Pfam-A with pyHMMER (gathering cutoffs by default), and
  resolve overlapping hits with cath-resolve-hits
- precomputed (import_annot_cds, enabled by --gff + --domtblout): reuse an
  existing genome-wide hmmsearch domtblout. CDSs falling inside each BGC are
  selected from the (tabix-indexed) GFF, the domtblout is subset per BGC, and
  the subsets are resolved with cath-resolve-hits
"""

###############################################################################
# 1. Set env
###############################################################################

import argparse
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from utils import (check_file, check_tools, extract_cds, find_files,
                   prepare_output_dir, resolve_domain_structure)

PFAM_DEFAULT = os.path.join(os.path.expanduser("~"), ".sedna", "db", "pfam", "Pfam-A.hmm")

###############################################################################
# 2. Define functions
###############################################################################

###############################################################################
# 2.1 Parse command-line arguments
###############################################################################

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"{os.path.basename(__file__)}: Protein-domain annotation of BGCs",
        add_help=False)

    parser.add_argument("--help", action="help", help="print this help")

    parser.add_argument("--input_dir", dest="input_dir", required=True,
        help="module 4 dereplicated/ directory (GBKs under gbks/, plus annot_metadata.tsv)")

    parser.add_argument("--pattern", dest="pattern", default=".*.gbk",
        help="regex for GBK file names (default: .*.gbk)")

    parser.add_argument("--cds_id", dest="cds_id", default="locus_tag",
        help="GBK qualifier used as CDS ID (default: locus_tag; protein_id for MIBiG-style records)")

    parser.add_argument("--hmm_db", dest="hmm_db", default=PFAM_DEFAULT,
        help=f"HMM database (default: {PFAM_DEFAULT})")

    parser.add_argument("--evalue_thres", dest="evalue_thres", type=float, default=1e-3,
        help="hmmsearch domain e-value threshold (default: 1e-3)")

    parser.add_argument("--cut_ga", dest="cut_ga", default=True,
        action=argparse.BooleanOptionalAction,
        help="use Pfam gathering cutoffs (default: --cut_ga)")

    parser.add_argument("--gff", dest="gff", default=None,
        help="precomputed mode: GFF of gene features of the sample's assembly")

    parser.add_argument("--domtblout", dest="domtblout", default=None,
        help="precomputed mode: genome-wide hmmsearch domtblout of the sample's proteins")

    parser.add_argument("--nslots", dest="nslots", type=int, default=4,
        help="number of threads (default: 4)")

    parser.add_argument("--output_dir", dest="output_dir", required=True,
        help="directory to output generated data")

    parser.add_argument("--overwrite", dest="overwrite", action="store_true", default=False,
        help="overwrite previous folder if present (default: False)")

    return parser.parse_args()

###############################################################################
# 2.2 annot_cds: extract CDSs of each GBK into FASTA files
###############################################################################

def create_fasta_files(paths2gbk: list, output_dir: str, cds_id: str) -> list:
    """Write <output_dir>/<gbk basename>.fasta per GBK; return the base paths."""
    path2names = []
    for gbk in paths2gbk:
        file_name = os.path.basename(os.path.splitext(gbk)[0])
        file_path = os.path.join(output_dir, file_name)
        extract_cds(input_gbk=gbk, cds_id=cds_id, output_fasta=f'{file_path}.fasta')
        path2names.append(file_path)
    return path2names

###############################################################################
# 2.3 annot_cds: hmmsearch with pyHMMER
###############################################################################

def run_hmm_annotation(path2names: list, hmm_db: str, evalue_thres: float,
                       nslots: int, cut_ga: bool) -> None:
    """Write <file>.domtblout for every <file>.fasta."""
    import pyhmmer

    bit_cutoffs = "gathering" if cut_ga else None
    for file in path2names:
        with pyhmmer.plan7.HMMFile(hmm_db) as hmms:
            with pyhmmer.easel.SequenceFile(f'{file}.fasta', digital=True) as seqs:
                with open(f'{file}.domtblout', "wb") as output_file:
                    for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, domE=evalue_thres,
                                                        bit_cutoffs=bit_cutoffs, cpus=nslots):
                        hits.write(output_file, format="domains", header=False)

###############################################################################
# 2.4 precomputed: sort, bgzip and tabix-index the GFF
###############################################################################

def sort_and_index_gff(gff: str, work_dir: str) -> str:
    """Return the path of a sorted, bgzipped, tabix-indexed copy of gff."""
    gff_name = os.path.splitext(os.path.basename(gff))[0]
    sorted_gff = os.path.join(work_dir, f'{gff_name}.sorted.gff.gz')

    if os.path.exists(f'{sorted_gff}.tbi'):
        return sorted_gff

    sort_gff_command = ["bash", "-c",
                        f"set -o pipefail; (grep '^#' {gff}; grep -v '^#' {gff} | "
                        f"sort -t \"$(printf '\\t')\" -k1,1 -k4,4n) | bgzip > {sorted_gff}"]
    result = subprocess.run(sort_gff_command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    # grep exits 1 when the GFF has no header lines; only a missing output is fatal
    if not os.path.isfile(sorted_gff) or os.path.getsize(sorted_gff) == 0:
        print(f"Error sorting GFF file:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    result = subprocess.run(["tabix", "-p", "gff", sorted_gff], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode != 0:
        print(f"Error indexing GFF file:\n{result.stderr}", file=sys.stderr)
        sys.exit(1)

    return sorted_gff

###############################################################################
# 2.5 precomputed: get the CDS IDs falling inside each BGC
###############################################################################

def extract_cds_ids(gff: str, annot_metadata_df) -> dict:
    """Return {bgc_id: [cds_id, ...]} from a tabix-indexed GFF."""
    import pysam

    cds_ids_dict = {}
    tabix_file = pysam.TabixFile(gff)

    for row in annot_metadata_df.itertuples(index=False):
        bgc_id = row.bgc_id
        start = row.start + 1  # convert coordinates to 1-based (as used in prodigal)
        region = f'{row.contig_id}:{start}-{row.end}'

        for record in tabix_file.fetch(region=region):
            record_list = record.split('\t')
            if record_list[2] == "CDS":
                cds_id = record_list[8].split(';')[0].replace('ID=', '')
                cds_ids_dict.setdefault(bgc_id, []).append(cds_id)

    tabix_file.close()
    return cds_ids_dict

###############################################################################
# 2.6 precomputed: subset the domtblout per BGC
###############################################################################

def filter_domtblout_by_cds(domtblout: str, cds_ids_dict: dict, output_dir: str) -> list:
    """Write <output_dir>/<bgc_id>.domtblout per BGC; return the base paths."""
    path2names = []
    for bgc_id, cds_ids in cds_ids_dict.items():
        cds_ids = set(cds_ids)
        base_path = os.path.join(output_dir, bgc_id)
        path2names.append(base_path)
        with open(domtblout, "r") as infile, open(f'{base_path}.domtblout', "w") as outfile:
            for line in infile:
                if not line.startswith("#"):
                    if line.strip().split()[0] in cds_ids:
                        outfile.write(line)
    return path2names

###############################################################################
# 3. Main
###############################################################################

def main() -> None:

    args = parse_args()
    precomputed = args.gff is not None or args.domtblout is not None

    ###########################################################################
    # 3.1 Sanity checks
    ###########################################################################

    check_tools(["cath-resolve-hits"])
    if not os.path.isdir(args.input_dir):
        print(f"Error: input dir {args.input_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    if precomputed:
        if args.gff is None or args.domtblout is None:
            print("Error: --gff and --domtblout must be given together", file=sys.stderr)
            sys.exit(1)
        check_tools(["bgzip", "tabix"])
        check_file(args.gff, f"--gff {args.gff}")
        check_file(args.domtblout, f"--domtblout {args.domtblout}")
    else:
        check_file(args.hmm_db, f"--hmm_db {args.hmm_db} (run bin/download_dbs.py first)")

    ###########################################################################
    # 3.2 Create output dir
    ###########################################################################

    prepare_output_dir(args.output_dir, args.overwrite)

    ###########################################################################
    # 3.3 Precomputed mode: import annotations
    ###########################################################################

    if precomputed:
        import pandas as pd

        annot_metadata = os.path.join(args.input_dir, "annot_metadata.tsv")
        if not os.path.isfile(annot_metadata):
            print(f"No annot_metadata.tsv in {args.input_dir}: no BGCs to annotate")
            return

        annot_metadata_df = pd.read_csv(annot_metadata, sep="\t")
        for col in ["bgc_id", "contig_id", "start", "end"]:
            if col not in annot_metadata_df.columns:
                print(f"Error: column '{col}' not found in {annot_metadata}", file=sys.stderr)
                sys.exit(1)

        sorted_gff = sort_and_index_gff(args.gff, os.getcwd())
        cds_ids_dict = extract_cds_ids(sorted_gff, annot_metadata_df)
        path2names = filter_domtblout_by_cds(args.domtblout, cds_ids_dict, args.output_dir)
        resolve_domain_structure(path2names)

        print(f"{os.path.basename(__file__)} executed successfully")
        return

    ###########################################################################
    # 3.4 Default mode: extract CDSs, hmmsearch, resolve
    ###########################################################################

    paths2gbk = find_files(input_dir=args.input_dir, pattern=args.pattern)
    if len(paths2gbk) == 0:
        print(f"No GBK files found in {args.input_dir}: no BGCs to annotate")
        return

    path2names = create_fasta_files(paths2gbk=paths2gbk, output_dir=args.output_dir,
                                    cds_id=args.cds_id)

    run_hmm_annotation(path2names=path2names, hmm_db=args.hmm_db,
                       evalue_thres=args.evalue_thres, nslots=args.nslots,
                       cut_ga=args.cut_ga)

    resolve_domain_structure(path2names)

    print(f"{os.path.basename(__file__)} executed successfully")


if __name__ == "__main__":
    main()
