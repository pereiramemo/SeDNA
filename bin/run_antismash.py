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

parser = argparse.ArgumentParser(prog='run_antismash.py', \
                                 description='Annotates BGC sequences utilizing the antiSMASH tool.')

# general parameters
parser.add_argument("--input_sample", required=True, help="Input fasta file.")
parser.add_argument("--threads", default = 4, type=int, help="Number of threads.")
parser.add_argument("--sample_name", default = "sample", help="Sample name.")
parser.add_argument("--output_dir", default = "sedna_output", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")

# antiSMASH parameters
parser.add_argument("--taxon", default = "bacteria", help="antiSMASH - {bacteria,fungi} Taxonomic classification of input sequence.")
parser.add_argument("--genefinding_tool", default = "prodigal-m", help="antiSMASH {glimmerhmm,prodigal,prodigal-m,none,error} Specify algorithm used for gene finding: GlimmerHMM, Prodigal, Prodigal Metagenomic/Anonymous mode, or none. The 'error' option will raise an error if genefinding is attempted. The 'none' option will not run genefinding.")
parser.add_argument("--minlength", type=int, default = 1000, help="antiSMASH - Only process sequences larger than <minlength>.")

# Get general parameters
args = parser.parse_args()
input_sample = args.input_sample
threads = int(args.threads)
output_dir = args.output_dir
sample_name = args.sample_name

# Get antiSMASH parameters
taxon = args.taxon
genefinding_tool = args.genefinding_tool
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

valid_taxa = {"bacteria", "fungi"}
if args.taxon not in valid_taxa:
    print(f"Error: Invalid taxon '{args.taxon}'. Must be one of {valid_taxa}.")
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
                
def run_antismash(input_sample, output_dir):
    command_antismash = [
                        "antismash", 
                         "--cpus", threads,
                         "--genefinding-tool", genefinding_tool,
                         "--taxon", taxon,
                         "--allow-long-headers",
                         "--minlength", minlength,
                         "--minimal",
                         "--output-dir", output_dir,
                         input_sample
                        ]

    result_antismash = subprocess.run(command_antismash, 
                                      shell=True, 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE, 
                                      text=True)

    if result_antismash.returncode != 0:
        print("Error executing antiSMASH")
        print("Error message:\n", result_antismash.stderr)
        sys.exit(1)
        
###############################################################################
#### 5 Create output dirs
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
#### 5.3. antismash output dir
###############################################################################

bgc_annot_inter_antismash_dir = os.path.join(bgc_annot_inter_dir, 'antismash')
create_directory(bgc_annot_inter_antismash_dir)

bgc_annot_inter_antismash_output_dir = os.path.join(bgc_annot_inter_antismash_dir, 'output')
create_directory(bgc_annot_inter_antismash_output_dir)

###############################################################################
## 6. Check if antismash output exists 
###############################################################################

bgc_annot_inter_antismash_dir_gbks = os.path.join(bgc_annot_inter_antismash_dir, 'gbks')

if os.path.exists(bgc_annot_inter_antismash_dir_gbks):
    print(f"Directory '{bgc_annot_inter_antismash_dir_gbks}' already exists. Use --overwrite to overwrite.")
    sys.exit()

###############################################################################
## 7. Run BGC annotation with antismash
###############################################################################

run_antismash(input_sample = input_sample, 
              output_dir = bgc_annot_inter_antismash_output_dir)

###############################################################################
## 8. Move files to the gbk folder
###############################################################################

# find all antismash GBK output files 
antismash_output_gbks = utilities.find_files(input_dir = bgc_annot_inter_antismash_output_dir,
                                             pattern = ".*.region.*.gbk")

if not antismash_output_gbks:
    print("No GBK files were generated by antiSMASH")
    print("No metadata file was generated after running antiSMASH")
    print("run_antismash.py executed successfully")
    sys.exit()

create_directory(bgc_annot_inter_antismash_dir_gbks)

for gbk in antismash_output_gbks:
        shutil.move(gbk, bgc_annot_inter_antismash_dir_gbks)

###############################################################################
## 9. Get metadata
###############################################################################

antismash_annot_metadata = utilities.antismash_annot_parser(input_dir = f'{bgc_annot_inter_antismash_dir}/gbks',
                                                            sample_name = sample_name,
                                                            input_fasta = input_sample)

###############################################################################
## 10. Export metadata
###############################################################################
                 
if antismash_annot_metadata is not None:
    antismash_annot_metadata_tsv = os.path.join(bgc_annot_inter_antismash_dir, "annot_metadata.tsv")
    antismash_annot_metadata.to_csv(antismash_annot_metadata_tsv, sep='\t', index=False)
else:
    print("No metadata file was generated after running antiSMASH")

###############################################################################
## 11. Print output message
###############################################################################

print("run_antismash.py executed successfully")