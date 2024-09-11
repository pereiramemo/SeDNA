#!/usr/bin/env python3

###############################################################################
## 1. Set env
###############################################################################

import argparse
import subprocess
import os
import sys
import sqlite3
import re
import shutil
from Bio import SeqIO
import import_ipynb
from src.utilities import *

###############################################################################
## 2. Parse input data and optional arguments
###############################################################################

parser = argparse.ArgumentParser(prog='bioprospecting.py', \
                                 description='Bioprospecting metagenomic data: BGCs identification, clustering and novelty assessment')

# general parameters
parser.add_argument("--input_sample", help="Input fasta file.")
parser.add_argument("--threads", default = 4, help="Number of threads.")
parser.add_argument("--sample_name", help="Sample name.")
parser.add_argument("--output_dir", help="Output directory.")
parser.add_argument("--overwrite", action="store_true", help="Ovewrite output directory.")
parser.add_argument("--antismash", action="store_true", help="Annotate BGCs with antiSMASH")
parser.add_argument("--deepbgc", action="store_true", help="Annotate BGCs with deepBGC")
parser.add_argument("--gecco", action="store_true", help="Annotate BGCs with gecco")

parser.add_argument("--antismash_output", default = None, help="Folder with precomputed anitSMASH BGC annotations")
parser.add_argument("--deepbgc_output",  default = None, help="Folder with precomputed deepBGC BGC annotations")
parser.add_argument("--gecco_output",  default = None, help="Folder with precomputed gecco BGC annotations")


# antiSMASH parameters
parser.add_argument("--taxon", default = "bacteria", help="antiSMASH - {bacteria,fungi} Taxonomic classification of input sequence.")
parser.add_argument("--genefinding_tool", default = "prodigal-m", help="antiSMASH {glimmerhmm,prodigal,prodigal-m,none,error} Specify algorithm used for gene finding: GlimmerHMM, Prodigal, Prodigal Metagenomic/Anonymous mode, or none. The 'error' option will raise an error if genefinding is attempted. The 'none' option will not run genefinding.")
parser.add_argument("--minlength", default = 1000, help="antiSMASH - Only process sequences larger than <minlength>")

# deepBGC parameters
parser.add_argument("--deepbgc_score_thres", default = 0.75, help = "Threshold value to filter out deepBGC annotated BGC sequences.")
parser.add_argument("--deepbgc_cds_count_thres", default = 2, help = "Threshold number of CDS to filter out deepBGC annotated BGC sequences.")
                    
# MMSeqs taxonomy parameters
parser.add_argument("--tax_lineage", default = 1, help="MMSeqs taxonomy - 0: don't show, 1: add all lineage names, 2: add all lineage taxids [0].")

# Get general paramteres
args = parser.parse_args()
input_sample = args.input_sample
threads = args.threads
output_dir = args.output_dir
input_bam = args.input_bam
sample_name = args.sample_name

# Get antiSMASH parameters
taxon = args.taxon
genefinding_tool = args.genefinding_tool
minlength = args.minlength
antismash_output = args.antismash_output

# Get deepBGC parameters
deepbgc_score_thres = float(args.deepbgc_score_thres)
deepbgc_cds_count_thres = int(args.deepbgc_cds_count_thres)

# Get deepBGC parameters
deepbgc_output = args.deepbgc_output

# Get gecco parameters
gecco_output = args.gecco_output

###############################################################################
## 3. Sanity checks and initializations
###############################################################################

# set annotation flags
antismash_annot_flag = False
deepbgc_annot_flag = False
gecco_annot_flag = False

# set metadata file names
antismash_annot_metadata_tsv = None
deepbgc_annot_metadata_tsv = None
gecco_annot_metadata_tsv = None

if args.antismash == False and args.deepbgc == False and args.gecco == False:
    if antismash_output is None and deepbgc_output is None and gecco_output is None:
        print(f"One BGC annotation tool must be selected or precomputed annotations provided")
        sys.exit()

###############################################################################
## 4. Create main output dir
###############################################################################

process = subprocess.Popen(
    "conda run -n ${CONDA_ENV_NAME} python script.py".split(), , stdout=subprocess.PIPE
)
output, error = process.communicate()

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
except FileExistsError:
    print(f"Directory '{output_dir}' already exists.")
    sys.exit()
except Exception as e:
    print(f'Error: {e}')
    sys.exit(1)

###############################################################################
#### 4.2 bgc_annot/inter/ output dir
###############################################################################

bgc_annot_inter_dir = f"{output_dir}/bgc_annot/inter"

try:
  os.makedirs(bgc_annot_inter_dir)
except Exception as e:
  print(f"An error occurred: {e}")
  sys.exit(1)

###############################################################################
#### 4.3 bgc_annot/sorted/ output dir
###############################################################################

bgc_annot_sorted_dir = f"{output_dir}/bgc_annot/sorted"

try:
  os.makedirs(bgc_annot_sorted_dir)
except Exception as e:
  print(f"An error occurred: {e}")
  sys.exit(1)

###############################################################################
#### 4.4 bgc_taxa output dir 
###############################################################################

bgc_taxa_dir = f"{output_dir}/bgc_taxa"

try:
  os.makedirs(bgc_taxa_dir)
except Exception as e:
  print(f"An error occurred: {e}")
  sys.exit(1)

###############################################################################
#### 4.5 output_tables dir
###############################################################################

output_tables_dir = f"{output_dir}/tables"

try:
  os.makedirs(output_tables_dir)
except Exception as e:
  print(f"An error occurred: {e}")
  sys.exit(1)

###############################################################################
## 5. BGC annotation
###############################################################################

###############################################################################
#### 5.1. Run antiSMASH
###############################################################################

if args.antismash == True and antismash_output is None:

    # create inter output
    bgc_annot_inter_antismash_dir = f"{bgc_annot_inter_dir}/antismash"
    try:
      os.makedirs(bgc_annot_inter_antismash_dir)
    except Exception as e:
      print(f"An error occurred: {e}")
      sys.exit(1)

    # # create sorted output
    # bgc_annot_sorted_antismash_dir = f"{bgc_annot_sorted_dir}/antismash"
    # try:
    #   os.makedirs(bgc_annot_sorted_antismash_dir)
    # except Exception as e:
    #   print(f"An error occurred: {e}")
    #   sys.exit(1)

    # run command
    current_directory = os.getcwd()
    antismash_output_current_dir = bgc_annot_inter_antismash_dir

    command_antismash = f"{current_directory}/src/run_antismash.sh \
                        {input_sample} {antismash_output_current_dir}  \
                          --cpus {threads} \
                          --genefinding-tool {genefinding_tool} \
                          --taxon {taxon} \
                          --allow-long-headers \
                          --minlength {minlength} \
                          --minimal"

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

    # set annotation flag: antismash was run 
    antismash_annot_flag = True


###############################################################################
#### 5.2. Run deepBGC
###############################################################################

if args.deepbgc == True and deepbgc_output is None:

    # create inter output
    bgc_annot_inter_deepbgc_dir = f"{bgc_annot_inter_dir}/deepbgc"
    try:
      os.makedirs(bgc_annot_inter_deepbgc_dir)
    except Exception as e:
      print(f"An error occurred: {e}")
      sys.exit(1)

    # # create sorted output
    # bgc_annot_sorted_deepbgc_dir = f"{bgc_annot_sorted_dir}/deepbgc"
    # try:
    #   os.makedirs(bgc_annot_sorted_deepbgc_dir)
    # except Exception as e:
    #   print(f"An error occurred: {e}")
    #   sys.exit(1)

    # run command
    current_directory = os.getcwd()
    deepbgc_output_current_dir = bgc_annot_inter_deepbgc_dir

    command_deepbgc = f"{current_directory}/src/run_deepbgc.sh \
                      {input_sample} \
                      {deepbgc_output_current_dir} \
                      {threads}"

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

    # set annotation flag: deepbgc was run 
    deepbgc_annot_flag = True

###############################################################################
#### 5.3. Run gecco
###############################################################################

if args.gecco == True and gecco_output is None:

    # create inter output
    bgc_annot_inter_gecco_dir = f"{bgc_annot_inter_dir}/gecco"
    try:
      os.makedirs(bgc_annot_inter_gecco_dir)
    except Exception as e:
      print(f"An error occurred: {e}")
      sys.exit(1)

    # # create sorted output
    # bgc_annot_sorted_gecco_dir = f"{bgc_annot_sorted_dir}/gecco"
    # try: 
    #   os.makedirs(bgc_annot_sorted_gecco_dir)
    # except Exception as e:
    #   print(f"An error occurred: {e}")
    #   sys.exit(1)

    # run command
    current_directory = os.getcwd()
    gecco_output_current_dir = bgc_annot_inter_gecco_dir

    command_gecco = f"{current_directory}/src/run_gecco.sh \
                    {input_sample} \
                    {gecco_output_current_dir}"

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

    # set annotation flag: gecco was run 
    gecco_annot_flag = True
