#! /usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import pandas as pd
import subprocess
import tempfile
import argparse
import os
import sys
import shutil

################################################################################
# 2. Parse parameters
################################################################################

parser = argparse.ArgumentParser(
    prog="get_bgc_coverage.py",
    description="Gets the coverage of BGC sequences from the MAGs associated BAM files."
)
parser.add_argument(
    "--annot_metadata", required=True,
    help="TSV file with BGC coordinates."
)

parser.add_argument(
    "--coverage", required = False, default = None,
    help="TSV file containing the coverage estimates per contig (generated with VEBA)."
)

parser.add_argument(
    "--bam_file", required = False, default = None, 
    help="BAM file containing the MAG coverage estimates."
)

parser.add_argument(
    "--overwrite", default = False, 
    action=argparse.BooleanOptionalAction, help="Overwrite output directory."
)

parser.add_argument(
    "--output_tsv", required = True,
    help="Output tsv file containing the BGC coverage information."
)
args = parser.parse_args()

annot_metadata = args.annot_metadata
bam_file = args.bam_file
coverage = args.coverage
overwrite = args.overwrite
output_tsv = args.output_tsv

################################################################################
# 3. Check if annot_metadata exists
################################################################################

if not os.path.exists(annot_metadata):
    print("annot_metadata file does not exist")
    sys.exit(0)
        
################################################################################
# 4. Check if the coverage or the bam file are provided
################################################################################

if coverage is None and bam_file is None:
    print("A coverage table or bam file must be provided")
    sys.exit(0)

if coverage is not None and bam_file is not None:
    print("If the coverage table is provided, the bam file will not be considered")
    sys.exit(0)

################################################################################
# 5. Check if output_tsv already exists
################################################################################

if os.path.exists(output_tsv):

    if args.overwrite is False:
        print(f'Output tsv {output_tsv} already exists. Use --overwrite to overwrite')
        sys.exit(0)

    if args.overwrite is True:
        try:
            os.remove(output_tsv)
        except OSError as e:
            print(f'Error: {e}')

################################################################################
# 6. Define functions
################################################################################

def parse_coverage_tables(coverage_df: pd.DataFrame, annot_metadata_df: pd.DataFrame) -> pd.DataFrame:
    """
    Filters rows in `coverage_df` where 'contigName' matches any 'contig_id' in `annot_metadata_df`,
    and includes the corresponding 'bgc_ids' from `annot_metadata_df`.

    Parameters:
    - coverage_df (pd.DataFrame): DataFrame with a 'contigName' column.
    - annot_metadata_df (pd.DataFrame): DataFrame with 'contig_id' and 'bgc_id' columns.

    Returns:
    - pd.DataFrame: Filtered DataFrame with an additional 'bgc_id' column.
    """
    
    # First, filter coverage_df using contig IDs
    contig_ids = annot_metadata_df['contig_id'].unique()
    filtered_df = coverage_df[coverage_df['contigName'].isin(contig_ids)]

    # Merge to bring in bgc_ids
    merged_df = pd.merge(
        filtered_df,
        annot_metadata_df[['contig_id', 'start', 'end', 'bgc_id']],
        how='left',
        left_on='contigName',
        right_on='contig_id'
    )

    # Drop unused columns
    merged_df = merged_df.drop(columns=['contig_id','contigLen','mapped.sorted.bam','mapped.sorted.bam-var'])
    
    # Re oder columns
    merged_df = merged_df[['contigName', 'start', 'end', 'bgc_id', 'totalAvgDepth']]

    return merged_df

def run_bedtools_coverage(bed_file, bam_file, output_tsv):
    """
    Runs the bedtools coverage command and writes output to a file.

    :param bed_file: Path to the BED file (regions of interest).
    :param bam_file: Path to the BAM file.
    :param output_file: Path to save the coverage output (default: coverage_summary.tsv).
    :return: The path of the output file.
    """
    try:
        # Construct the command
        command = f"bedtools coverage -a {bed_file} -b {bam_file} -mean > {output_tsv}"

        # Run the command
        subprocess.run(command, shell=True, check=True)
        return output_tsv

    except subprocess.CalledProcessError as e:
        print(f"Error running bedtools coverage: {e}")
        sys.exit(1)

################################################################################
# 7. Import annot_metadata table
################################################################################

annot_metadata_df = pd.read_csv(annot_metadata, sep="\t")

################################################################################
# 8. Import coverage table, parse coverage estimates and write output
################################################################################

if coverage is not None:
    
    coverage_df = pd.read_csv(coverage, sep="\t")
    output_df = parse_coverage_tables(coverage_df = coverage_df,
                                       annot_metadata_df = annot_metadata_df)
    output_df.to_csv(output_tsv, sep="\t", index=False, header = False)

################################################################################
# 9. Create and write BED file and estimate coverage
################################################################################
 
if bam_file is not None and coverage is None:

    bed_df = annot_metadata_df[['contig_id', 'start', 'end', 'bgc_id']]
    input_dir = os.path.dirname(annot_metadata)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".bed", mode='w') as tmp_bed:
        bed_df.to_csv(tmp_bed.name, sep="\t", index=False, header=False)
        bed_file = tmp_bed.name

    run_bedtools_coverage(bed_file = bed_file, 
                          bam_file = bam_file, 
                          output_tsv = output_tsv)

################################################################################
# 10. End of execution
################################################################################

print("get_bgc_coverage module executed successfully")
