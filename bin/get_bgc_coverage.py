#! /usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import pandas as pd
import subprocess
import argparse
import os
import sys
import shutil
import re

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
    "--bam_file", required=True,
    help="BAM file containing the MAG coverage estimates."
)

parser.add_argument(
    "--overwrite", default = False, 
    action=argparse.BooleanOptionalAction, help="Overwrite output directory."
)

parser.add_argument(
    "--output_tsv",
    help="Output tsv file containing the BGC coverage information."
)
args = parser.parse_args()

bam_file = args.bam_file
annot_metadata = args.annot_metadata
overwrite = args.overwrite
output_tsv = args.output_tsv

################################################################################
# 3. Check if annot_metadata exists
################################################################################

if not os.path.exists(annot_metadata):
    print("annot_metadata file does not exists")
    sys.exit(0)

################################################################################
# 4. Check if output_tsv already exists
################################################################################

if os.path.exists(output_tsv) is True:

    if args.overwrite is False:
        print(f'Output tsv {output_tsv} already exists. Use --overwrite to overwrite')
        sys.exit(0)

    if args.overwrite is True:
        try:
            os.remove(output_tsv)
        except OSError as e:
            print(f'Error: {e}')

################################################################################
# 5. Define functions
################################################################################

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
# 6. Import annot_metadata.tsv
################################################################################

annot_metadata_df = pd.read_csv(annot_metadata, sep="\t")

################################################################################
# 7. Check if annot_metadata already has coverage estimates
################################################################################

if 'coverage' in annot_metadata_df.columns:
  print(f'{annot_metadata} already includes the coverage estimate')
  sys.exit(0)

################################################################################
# 8. Create BED file
################################################################################

bed_df = annot_metadata_df[['contig_id', 'start', 'end', 'bgc_id']]

################################################################################
# 9. Create and write BED file
################################################################################

input_dir = os.path.dirname(annot_metadata)
bed_file = os.path.join(input_dir,"bgcs.bed")
bed_df.to_csv(bed_file, sep = "\t",  index=False, header=False)

################################################################################
# 10. Estimate coverage 
################################################################################

run_bedtools_coverage(bed_file = bed_file, 
                      bam_file = bam_file, 
                      output_tsv = output_tsv)

################################################################################
# 11. Add coverage to annot_metadata 
################################################################################

bgc_coverage_df = pd.read_csv(output_tsv, sep="\t", header = None)
bgc_coverage_df.columns = ['contig_id', 'start', 'end', 'bgc_id', "coverage"]

annot_metadata_ext_df = pd.merge(annot_metadata_df, bgc_coverage_df,
                                 on=['contig_id', 'start', 'end', 'bgc_id'], 
                                 how='left')

################################################################################
# 12. Rename annot_metadata  
################################################################################

input_dir = os.path.dirname(annot_metadata)
annot_metadata_bkup = os.path.join(input_dir, "annot_metadata_bkup.tsv")
os.rename(annot_metadata, annot_metadata_bkup)

################################################################################
# 13. Export annot_metadata_ext  
################################################################################

input_dir = os.path.dirname(annot_metadata)
annot_metadata_ext_file = os.path.join(input_dir, "annot_metadata.tsv")
annot_metadata_ext_df.to_csv(annot_metadata_ext_file, sep = "\t", 
                             index=False, header=True)

################################################################################
# 14. Test that annot_metadata_ext_file exists and is not empty
################################################################################

if os.path.exists(annot_metadata_ext_file) and os.path.getsize(annot_metadata_ext_file) > 0:
  os.remove(annot_metadata_bkup)
else:
  print(f"File {annot_metadata_ext_file} not found or empty")
  os.rename(annot_metadata_bkup, annot_metadata)
  sys.exit(1)

################################################################################
# 15. End of execution
################################################################################

print("get_bgc_coverage module executed successfully")
