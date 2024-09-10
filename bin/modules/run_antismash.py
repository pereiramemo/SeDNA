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

parser = argparse.ArgumentParser(prog='run_antismash.py', \
                                 description='Annotates BGC sequences utilizing the antiSMASH tool')

# general parameters
parser.add_argument("--input_sample", help="Input fasta file.")
parser.add_argument("--threads", default = 4, help="Number of threads.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--output_dir", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Ovewrite output directory.")

# antiSMASH parameters
parser.add_argument("--taxon", default = "bacteria", help="antiSMASH - {bacteria,fungi} Taxonomic classification of input sequence.")
parser.add_argument("--genefinding_tool", default = "prodigal-m", help="antiSMASH {glimmerhmm,prodigal,prodigal-m,none,error} Specify algorithm used for gene finding: GlimmerHMM, Prodigal, Prodigal Metagenomic/Anonymous mode, or none. The 'error' option will raise an error if genefinding is attempted. The 'none' option will not run genefinding.")
parser.add_argument("--minlength", default = 1000, help="antiSMASH - Only process sequences larger than <minlength>")

# Get general paramteres
args = parser.parse_args()
input_sample = args.input_sample
threads = args.threads
output_dir = args.output_dir
sample_name = args.sample_name

# Get antiSMASH parameters
taxon = args.taxon
genefinding_tool = args.genefinding_tool
minlength = args.minlength

###############################################################################
## 3. Sanity checks and initializations
###############################################################################

# set metadata file names
antismash_annot_metadata_tsv = None

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
#### 4.1 bgc_annot/inter/ output dir
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
## 5.Create antismash output dir
###############################################################################

# create inter output
bgc_annot_inter_antismash_dir = f"{bgc_annot_inter_dir}/antismash"
try:
    os.makedirs(bgc_annot_inter_antismash_dir)
except Exception as e:
    print(f"Error: {e}")
    sys.exit(1)

###############################################################################
## 6. Run BGC annotation with antismash
###############################################################################

    
# # run command
current_directory = os.getcwd()
antismash_output_current_dir = bgc_annot_inter_antismash_dir

command_antismash = f"antismash \
                       --cpus {threads} \
                       --genefinding-tool {genefinding_tool} \
                       --taxon {taxon} \
                       --allow-long-headers \
                       --minlength {minlength} \
                       --minimal \
                       --output-dir {antismash_output_current_dir}  \
                        {input_sample}" 

result_antismash = subprocess.run(command_antismash, 
                                  shell=True, 
                                  stdout=subprocess.PIPE, 
                                  stderr=subprocess.PIPE, 
                                  text=True)

if result_antismash.returncode == 0:
    print("antiSMASH executed successfully")
else:
    print("Error executing antiSMASH")
    print("Error message:\n", result_antismash.stderr)
    sys.exit()

###############################################################################
#### 6 anitsmash BGC metadata
###############################################################################

# find all antismash GBK output files 
antismash_output_gbks = utilities.find_files(input_dir = antismash_output_current_dir,
                                             pattern = ".*.region.*.gbk")

try:
    os.makedirs(f'{antismash_output_current_dir}/gbks')
except FileExistsError:
    print(f"Directory '{antismash_output_current_dir}/gbks' already exists.")
    sys.exit()

for gbk in antismash_output_gbks:
        shutil.move(gbk, f'{antismash_output_current_dir}/gbks')

antismash_annot_metadata = utilities.antismash_annot_parser(input_dir = f'{antismash_output_current_dir}/gbks',
                                                            sample_name = sample_name,
                                                            input_fasta = input_sample)

if antismash_annot_metadata is not None:
    antismash_annot_metadata_tsv = f'{bgc_annot_inter_antismash_dir}/antismash_annot_metadata.tsv'
    antismash_annot_metadata.to_csv(antismash_annot_metadata_tsv, sep='\t', index=False)
else:
    antismash_annot_metadata_tsv = None
