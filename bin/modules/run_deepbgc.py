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

current_dir = os.getcwd()
sys.path.append(f'{current_dir}/src')
from bgc_annot import utilities

###############################################################################
## 2. Parse input data and optional arguments
###############################################################################

parser = argparse.ArgumentParser(prog='run_deepbgc.py', \
                                 description='Annotates BGC sequences utilizing the deepBGC tool')

# general parameters
parser.add_argument("--input_sample", help="Input fasta file.")
parser.add_argument("--threads", default = 4, help="Number of threads.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--output_dir", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Ovewrite output directory.")
# deepBGC parameters
parser.add_argument("--deepbgc_score_thres", default = 0.75, help = "Threshold value to filter out deepBGC annotated BGC sequences.")
parser.add_argument("--deepbgc_cds_count_thres", default = 2, help = "Threshold number of CDS to filter out deepBGC annotated BGC sequences.")


# Get general parameters
args = parser.parse_args()
input_sample = args.input_sample
threads = args.threads
output_dir = args.output_dir
sample_name = args.sample_name

# Get deepBGC parameters
deepbgc_score_thres = float(args.deepbgc_score_thres)
deepbgc_cds_count_thres = int(args.deepbgc_cds_count_thres)

###############################################################################
## 3. Sanity checks and initializations
###############################################################################

# set metadata file names
deepbgc_annot_metadata_tsv = None

###############################################################################
#### 4 Create output dirs
###############################################################################

###############################################################################
#### 4.1 Main output dir
###############################################################################

output_dir = f'{output_dir}/{sample_name}'

if args.overwrite and os.path.exists(output_dir):
    try:
        shutil.rmtree(output_dir)
    except OSError as e:
        print(f'Error: {e}')
    try:       
        os.makedirs(output_dir)
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)

###############################################################################
#### 4.2 bgc_annot/inter/ output dir
###############################################################################

bgc_annot_inter_dir = f"{output_dir}/bgc_annot/inter"

if not os.path.exists(output_dir):
    try:
        os.makedirs(bgc_annot_inter_dir)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

###############################################################################
#### 4.3 bgc_annot/sorted/ output dir
###############################################################################

bgc_annot_sorted_dir = f"{output_dir}/bgc_annot/sorted"

if not os.path.exists(bgc_annot_sorted_dir):
    try:
        os.makedirs(bgc_annot_sorted_dir)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

###############################################################################
## 4.4. deepbgc output dir
###############################################################################

# create inter output
bgc_annot_inter_deepbgc_dir = f"{bgc_annot_inter_dir}/deepbgc"
try:
    os.makedirs(bgc_annot_inter_deepbgc_dir)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
    
bgc_annot_inter_deepbgc_output_dir = f"{bgc_annot_inter_deepbgc_dir}/output"
try:
    os.makedirs(bgc_annot_inter_deepbgc_output_dir)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

###############################################################################
## 5. Run BGC annotation with deepbgc
###############################################################################

# run command
command_deepbgc = f"deepbgc pipeline \
                  --detector deepbgc \
                  --classifier product_class \
                  --prodigal-meta-mode \
                  --classifier-score {deepbgc_score_thres} \
                  --min-domains {deepbgc_cds_count_thres} \
                  --output {bgc_annot_inter_deepbgc_output_dir} \
                  {input_sample}" 

result_deepbgc = subprocess.run(command_deepbgc, 
                                  shell=True, 
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE, 
                                  text=True)

if result_deepbgc.returncode == 0:
    print("deepBGC executed successfully")
else:
    print("Error executing deepBGC")
    print("Error message:\n", result_deepbgc.stderr)
    sys.exit()

###############################################################################
#### 6 deepBGC BGC metadata
###############################################################################

# find the sinlge GBK output file from deepBGC
deepbgc_output_gbk = utilities.find_files(input_dir = bgc_annot_inter_deepbgc_output_dir, 
                                          pattern = ".bgc.gbk")

if len(deepbgc_output_gbk) > 1:
    print("Error message:\nMore than one deepBGC GBK file output found\nThere should be only one.")
    print(deepbgc_output_gbk)
    sys.exit() 

if len(deepbgc_output_gbk) == 1:
    # split multiple GBK into different files
    utilities.gbk_splitter(input_file = deepbgc_output_gbk[0],
                           output_dir = f'{bgc_annot_inter_deepbgc_dir}/gbks', 
                           sample_name = sample_name)

    # Filter deebBGC outputs with low score
    utilities.gbk_filter(input_dir = f'{bgc_annot_inter_deepbgc_dir}/gbks',
                         sample_name = sample_name,
                         deepbgc_score_thres = deepbgc_score_thres,
                         deepbgc_cds_count_thres = deepbgc_cds_count_thres,
                         output_dir = f'{bgc_annot_inter_deepbgc_dir}/gbks_removed')

# Get metadata
deepbgc_annot_metadata = utilities.deepbgc_annot_parser(input_dir = f'{bgc_annot_inter_deepbgc_dir}/gbks',
                                                        sample_name = sample_name,
                                                        input_fasta = input_sample)


# Export metadata
if deepbgc_annot_metadata is not None:
    deepbgc_annot_metadata_tsv = f'{bgc_annot_inter_deepbgc_dir}/deepbgc_annot_metadata.tsv'
    deepbgc_annot_metadata.to_csv(deepbgc_annot_metadata_tsv, sep='\t', index=False)
else:
    deepbgc_annot_metadata_tsv = None
