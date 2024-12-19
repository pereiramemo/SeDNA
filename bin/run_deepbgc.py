#!/usr/bin/env python3

###############################################################################
## 1. Set env
###############################################################################

import argparse
import subprocess
import os
import sys
import re
import shutil
from Bio import SeqIO

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{current_dir}/src')
from bgc_annot import utilities

###############################################################################
## 2. Parse input data and optional arguments
###############################################################################

parser = argparse.ArgumentParser(prog='run_deepbgc.py', \
                                 description='Annotates BGC sequences utilizing the deepBGC tool.')

# general parameters
parser.add_argument("--input_sample", required=True, help="Input fasta file.")
parser.add_argument("--threads", default = 4, type=int, help="Number of threads.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--output_dir", default = "sedna_output", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")
# deepBGC parameters
parser.add_argument("--score_thres", default = 0.80,  type=float, help = "deepBGC - Threshold value to filter out deepBGC annotated BGC sequences.")
parser.add_argument("--cds_count_thres", default = 3,  type=int, help = "deepBGC - Threshold number of CDS to filter out deepBGC annotated BGC sequences.")
parser.add_argument("--minlength", default = 1000,  type=int, help = "deepBGC - Minimum BGC nucleotide length")

# Get general parameters
args = parser.parse_args()
input_sample = args.input_sample
threads = int(args.threads)
output_dir = args.output_dir
sample_name = args.sample_name

# Get deepBGC parameters
score_thres = float(args.score_thres)
cds_count_thres = int(args.cds_count_thres)
minlength = int(args.minlength)

###############################################################################
## 3. Sanity checks
###############################################################################

try:
    with open(input_sample, "r") as handle:
        fasta_check = any(SeqIO.parse(handle, "fasta"))
    if not fasta_check:
        raise ValueError(f"{input_sample} contains no valid FASTA records.")
except (FileNotFoundError, ValueError) as e:
    print(f"Error: {e}")
    sys.exit(1)
    
if threads <= 0:
    print(f"Error: Invalid number of threads. Must be a positive integer.")
    sys.exit(1)
    
if minlength <= 0:
    print(f"Error: Invalid minimum length. Must be a positive integer.")
    sys.exit(1)    

###############################################################################
## 4. Define functions
###############################################################################

def remove_directory(path):
    if os.path.exists(path):
        try:
            shutil.rmtree(path)
        except OSError as e:
            print(f'Error: {e}')

def create_directory(path):
    try:
        os.makedirs(path, exist_ok=True)
    except Exception as e:
        print(f"Error creating directory {path}: {e}")
        sys.exit(1)

import subprocess

def run_deepbgc(input_sample, output_dir, minlength = minlength):
    """
    Run the DeepBGC pipeline with the specified parameters.

    Args:
        input_sample (str): Path to the input sample file.
        output_dir (str): Directory to save the output.
        minlength (int): Minimum nucleotide length for processing.

    Raises:
        RuntimeError: If the DeepBGC pipeline execution fails.
    """
    # Ensure minlength is passed as a string
    command_deepbgc = [
        "deepbgc", "pipeline",
        "--detector", "deepbgc",
        "--classifier", "product_class",
        "--prodigal-meta-mode",
        "--min-nucl", str(minlength),
        "--output", output_dir,
        input_sample
    ]

    try:
        # Run the command
        result_deepbgc = subprocess.run(
            command_deepbgc,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True  # Raises CalledProcessError if the command fails
        )

    except subprocess.CalledProcessError as e:
        # Handle and raise a more descriptive error
        error_msg = (
            f"Error executing DeepBGC pipeline:\n"
            f"Command: {' '.join(command_deepbgc)}\n"
            f"Error Message: {e.stderr.strip()}"
        )
        raise RuntimeError(error_msg) from e

###############################################################################
## 5. Create output dirs
###############################################################################

###############################################################################
#### 5.1 Main output dir
###############################################################################

output_dir = os.path.join(output_dir, sample_name)

if os.path.exists(output_dir):
    if args.overwrite:
        remove_directory(output_dir)
        create_directory(output_dir)
else:
    create_directory(output_dir)

###############################################################################
#### 5.2 bgc_annot/inter/ output dir
###############################################################################

bgc_annot_inter_dir = os.path.join(output_dir, "bgc_annot", "inter")
create_directory(bgc_annot_inter_dir)

###############################################################################
#### 5.3. deepbgc output dir
###############################################################################

bgc_annot_inter_deepbgc_dir = os.path.join(bgc_annot_inter_dir, 'deepbgc')
create_directory(bgc_annot_inter_deepbgc_dir)

bgc_annot_inter_deepbgc_output_dir = os.path.join(bgc_annot_inter_deepbgc_dir, 'output')
create_directory(bgc_annot_inter_deepbgc_output_dir)

###############################################################################
## 6. Check if deepbgc output exists 
###############################################################################

bgc_annot_inter_deepbgc_dir_gbks = os.path.join(bgc_annot_inter_deepbgc_dir, 'gbks')

if os.path.exists(bgc_annot_inter_deepbgc_dir_gbks):
    print(f"Directory '{bgc_annot_inter_deepbgc_dir_gbks}' already exists. Use --overwrite to overwrite.")
    sys.exit()

###############################################################################
## 7. Run BGC annotation with deepbgc
###############################################################################

# run command
run_deepbgc(input_sample = input_sample, 
            output_dir = bgc_annot_inter_deepbgc_output_dir)

###############################################################################
## 8. Find GBK output (single) file from deepBGC
###############################################################################

deepbgc_output_gbk = utilities.find_files(input_dir = bgc_annot_inter_deepbgc_output_dir, 
                                          pattern = ".bgc.gbk")

###############################################################################
## 9. Split multiple GBK into different files
###############################################################################

if len(deepbgc_output_gbk) > 1:
    print("Error message:\nMore than one deepBGC GBK file output found\nThere should be only one.")
    print(deepbgc_output_gbk)
    sys.exit(1)

if len(deepbgc_output_gbk) == 0:
    print("No GBK files were generated by deepBGC")
    print("No metadata file was generated after running deepBGC")
    print("run_antismash.py executed successfully")
    sys.exit()
    
create_directory(bgc_annot_inter_deepbgc_dir_gbks)

utilities.gbk_splitter(input_file = deepbgc_output_gbk[0],
                       output_dir = bgc_annot_inter_deepbgc_dir_gbks, 
                       sample_name = sample_name)

###############################################################################
## 10. Filter deepBGC outputs with low score
###############################################################################

bgc_annot_inter_deepbgc_dir_gbks_removed = os.path.join(bgc_annot_inter_deepbgc_dir, 'gbks_removed')

utilities.gbk_filter(input_dir = bgc_annot_inter_deepbgc_dir_gbks,
                     sample_name = sample_name,
                     deepbgc_score_thres = score_thres,
                     deepbgc_cds_count_thres = cds_count_thres,
                     output_dir = bgc_annot_inter_deepbgc_dir_gbks_removed)

# Note: this filtering could have been done when running the tool; however, to be avoid running the tool
# more than once to reduce stringency, this is done separately and the full output is kept (as obtained with default values).

###############################################################################
## 11. Get metadata
###############################################################################

deepbgc_annot_metadata = utilities.deepbgc_annot_parser(input_dir = bgc_annot_inter_deepbgc_dir_gbks,
                                                        sample_name = sample_name,
                                                        input_fasta = input_sample)

###############################################################################
## 12. Export metadata
###############################################################################

if deepbgc_annot_metadata is not None:
    deepbgc_annot_metadata_tsv = os.path.join(bgc_annot_inter_deepbgc_dir,'annot_metadata.tsv')
    deepbgc_annot_metadata.to_csv(deepbgc_annot_metadata_tsv, sep='\t', index=False)
else:
    print("No metadata file was generated after running deepBGC")

###############################################################################
## 13. Format deepBGC GBKs
###############################################################################

if deepbgc_annot_metadata is not None:
    
    utilities.format_gbks(input_dir = bgc_annot_inter_deepbgc_dir_gbks,
                          input_tsv = deepbgc_annot_metadata_tsv,
                          output_dir = bgc_annot_inter_deepbgc_dir_gbks,
                          tool = "deepbgc")

###############################################################################
## 14. Print output message
###############################################################################
    
print("run_deepbgc.py executed successfully")