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

parser = argparse.ArgumentParser(prog='run_gecco.py', \
                                 description='Annotates BGC sequences utilizing the gecco tool.')

# general parameters
parser.add_argument("--input_sample", help="Input fasta file.")
parser.add_argument("--threads", default = 4, type = int, help="Number of threads.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--output_dir",  default = "sedna_output", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")

# Get general parameters
args = parser.parse_args()
input_sample = args.input_sample
threads = int(args.threads)
output_dir = args.output_dir
sample_name = args.sample_name

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

# if minlength <= 0:
#     print(f"Error: Invalid minimum length. Must be a positive integer.")
#     sys.exit(1)    

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

def run_gecco(input_sample, output_dir):

    command_gecco = f"gecco run \
                    --genome {input_sample} \
                    --jobs {threads} \
                    --output-dir {output_dir}"

    result_gecco = subprocess.run(command_gecco, 
                                  shell=True, 
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE, 
                                  text=True)

    if result_gecco.returncode != 0:
        print("Error executing gecco")
        print("Error message:\n", result_gecco.stderr)
        sys.exit()


def move_gbks(gecco_gbk_list, output_dir):
    for file_path in gecco_gbk_list:
        file = os.path.basename(file_path)
        pattern = r'_cluster_(\d+)\.gbk'
        match = re.search(pattern, file)
        number = int(match.group(1))
        file_renamed = re.sub(r'_cluster_\d+.gbk','.region{:03d}.gbk'.format(number), file)
        file_path_renamed = os.path.join(output_dir, file_renamed)
        shutil.copy(file_path, file_path_renamed)

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
#### 5.3 gecco output dir
###############################################################################

bgc_annot_inter_gecco_dir = os.path.join(bgc_annot_inter_dir, 'gecco')
create_directory(bgc_annot_inter_gecco_dir)

bgc_annot_inter_gecco_output_dir = os.path.join(bgc_annot_inter_gecco_dir, 'output')
create_directory(bgc_annot_inter_gecco_output_dir)

###############################################################################
## 6. Check if deepbgc output exists 
###############################################################################

bgc_annot_inter_gecco_dir_gbks = os.path.join(bgc_annot_inter_gecco_dir, 'gbks')

if os.path.exists(bgc_annot_inter_gecco_dir_gbks):
    print(f"Directory '{bgc_annot_inter_gecco_dir_gbks}' already exists. Use --overwrite to overwrite.")
    sys.exit()

###############################################################################
## 7. Run BGC annotation with gecco
###############################################################################

run_gecco(input_sample, bgc_annot_inter_gecco_output_dir)

###############################################################################
## 8. Move files to the gbk folder
###############################################################################

# copy all gecco output gbk files into a gbk folder and rename these following \
# the format region*.gbk (as antiSMASH).
gecco_gbk_list = utilities.find_files(input_dir =  bgc_annot_inter_gecco_output_dir, 
                                      pattern = r'.cluster_[0-9]+.gbk')

create_directory(bgc_annot_inter_gecco_dir_gbks)

move_gbks(gecco_gbk_list = gecco_gbk_list, 
          output_dir = bgc_annot_inter_gecco_dir_gbks)

###############################################################################
## 9. Get metadata
###############################################################################

# Note that to generate the metadata table for gecco, the output tsv table generated by
# gecco is used, and thus, there is no need to parse GBK files. Accordingly, the input
# directory is not the gbk folder.
gecco_annot_metadata = utilities.gecco_annot_parser(input_dir = bgc_annot_inter_gecco_output_dir,
                                                    sample_name = sample_name,
                                                    input_fasta = input_sample)

###############################################################################
## 10. Export metadata 
###############################################################################

if gecco_annot_metadata is not None:
    gecco_annot_metadata_tsv = f'{bgc_annot_inter_gecco_dir}/annot_metadata.tsv'
    gecco_annot_metadata.to_csv(gecco_annot_metadata_tsv, sep='\t', index=False)
else:
    print("No metadata file was generated after running GECCO")

###############################################################################
## 11. Print output message 
###############################################################################

print("run_gecco.py executed successfully")