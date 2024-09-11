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

parser = argparse.ArgumentParser(prog='run_gecco.py', \
                                 description='Annotates BGC sequences utilizing the gecco tool')

# general parameters
parser.add_argument("--input_sample", help="Input fasta file.")
parser.add_argument("--threads", default = 4, help="Number of threads.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--output_dir", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Ovewrite output directory.")

# Get general parameters
args = parser.parse_args()
input_sample = args.input_sample
threads = args.threads
output_dir = args.output_dir
sample_name = args.sample_name

###############################################################################
## 3. Sanity checks and initializations
###############################################################################

# set metadata file names
gecco_annot_metadata_tsv = None

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
## 4.4. gecco output dir
###############################################################################

# create inter output
bgc_annot_inter_gecco_dir = f"{bgc_annot_inter_dir}/gecco"
try:
    os.makedirs(bgc_annot_inter_gecco_dir)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)
    
bgc_annot_inter_gecco_output_dir = f"{bgc_annot_inter_gecco_dir}/output"
try:
    os.makedirs(bgc_annot_inter_gecco_output_dir)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)    

###############################################################################
## 5. Run BGC annotation with gecco
###############################################################################

# run command
command_gecco = f"gecco run \
                   --genome {input_sample} \
                   --jobs {threads} \
                   --output-dir {bgc_annot_inter_gecco_output_dir}"

result_gecco = subprocess.run(command_gecco, 
                              shell=True, 
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.PIPE, 
                              text=True)

if result_gecco.returncode == 0:
    print("gecco executed successfully")
else:
    print("Error executing gecco")
    print("Error message:\n", result_gecco.stderr)
    sys.exit()

###############################################################################
#### 6. gecco BGC metadata
###############################################################################

# copy all gecco output gbk files into a gbk folder and rename these following \
# the format region*.gbk (as antiSMASH).
gecco_gbk_list = utilities.find_files(input_dir =  bgc_annot_inter_gecco_output_dir, 
                                      pattern = r'.cluster_[0-9]+.gbk')

if not os.path.exists(f'{bgc_annot_inter_gecco_dir}/gbks'):
    os.makedirs(f'{bgc_annot_inter_gecco_dir}/gbks')


for file_path in gecco_gbk_list:
    file = os.path.basename(file_path)
    pattern = r'_cluster_(\d+)\.gbk'
    match = re.search(pattern, file)
    number = int(match.group(1))
    file_renamed = re.sub(r'_cluster_\d+.gbk','.region{:03d}.gbk'.format(number), file)
    file_path_renamed = os.path.join(f'{bgc_annot_inter_gecco_dir}/gbks', file_renamed)
    shutil.copy(file_path, file_path_renamed)
    
# Get metadata
# Note that to generate the metadata table for gecco, the output tsv table generated by
# this tool is used, and thus, there is no need to parse GBK files. Accordingly, the input
# directory is not the gbk folder.
gecco_annot_metadata = utilities.gecco_annot_parser(input_dir = bgc_annot_inter_gecco_output_dir,
                                                    sample_name = sample_name,
                                                    input_fasta = input_sample)
# Export metadata 
if gecco_annot_metadata is not None:
    gecco_annot_metadata_tsv = f'{bgc_annot_inter_gecco_dir}/annot_metadata.tsv'
    gecco_annot_metadata.to_csv(gecco_annot_metadata_tsv, sep='\t', index=False)
else:
    gecco_annot_metadata_tsv = None
