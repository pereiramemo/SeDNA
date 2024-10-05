#! /usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import pandas as pd
import os
import re
import csv
import shutil
import argparse
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir = os.path.dirname(current_dir)
sys.path.append(f'{bin_dir}/src')
from bgc_derep import dsc, dnsc, outputs 

################################################################################
# 2. Parse parameters
################################################################################

parser = argparse.ArgumentParser(prog='dereplicate_bgcs.py', 
                                 description='Dereplicate BGCs annotated with different tools.')
parser.add_argument("--input_dir", default = None, help = "The input directory having the metadata tables (tsv) generated with the run_* modules.")
parser.add_argument("--overlap_thres", default = "0.75", 
                    help="Percentage of overlap (in relation to longest BGCs) to determine if two BGCs are overlapped or partially overlapped.")
parser.add_argument("--metadata", default = None, help="Comma separated list of tsv tables containing the metadata of annotated BGCs.")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")
parser.add_argument("--output_dir", help = "The output directory where non-overlapped, partially overlapped and overlapped BGCs sequences and metadata tables (tsv).")

args = parser.parse_args()
input_dir = args.input_dir
overlap_thres = args.overlap_thres
metadata = args.metadata
overwrite = args.overwrite
output_dir = args.output_dir

################################################################################
# 3. Create or remove output directory
################################################################################

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
        
if args.overwrite is False and os.path.exists(output_dir) is True:
    print(f'Output dir {output_dir} already exists. Use --overwrite to overwrite')
    sys.exit(0)

################################################################################
# 4. Load the metadata tables and check if exist
################################################################################

if metadata is None and input_dir is not None:
    metadata_list = list()
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if re.search("annot_metadata.tsv", file):
                metadata_list.append(os.path.join(root, file))
    
else:
    metadata_list = [i.strip() for i in metadata.split(',')]
    
metadata_dict = dict()
for i in range(0, len(metadata_list)):
    file = metadata_list[i]
    if not os.path.exists(file):
        print(f'File {file} not found')
        sys.exit(1)
    else:
        metadata_df = pd.read_csv(file, sep='\t')
        tool = metadata_df['tool'].iloc[0]
        metadata_dict[tool] = metadata_df
        metadata_list[i] = metadata_df
    
################################################################################
# 5. Initialize dictionaries
################################################################################

dereplicated_bgcs = dict()
non_overlapped_bgcs = dict()
partially_overlapped_bgcs = dict()
overlapped_bgcs = dict()

for tool in metadata_dict:
    dereplicated_bgcs[tool] = dict()
    non_overlapped_bgcs[tool] = dict()
    partially_overlapped_bgcs[tool] = dict()
    overlapped_bgcs[tool] = dict()

###############################################################################
## 6. Define recursive function
###############################################################################

def recursive_dereplication(metadata: list =  None,
                            dereplicated_bgcs: str = None,
                            overlapped_bgcs: str = None,
                            partially_overlapped_bgcs: str = None,
                            non_overlapped_bgcs: str = None):

    metadata_range = list(range(0,len(metadata)))
    for i in metadata_range:

        current_i = i
        next_i = i +1
        metadata_df1 = metadata[current_i]
        metadata_df2 = metadata[next_i]

        output_shared_contigs = dsc.dereplicate_shared_contigs(metadata_df1= metadata_df1, 
                                                               metadata_df2 = metadata_df2,
                                                               dereplicated_bgcs = dereplicated_bgcs,
                                                               overlapped_bgcs = overlapped_bgcs,
                                                               partially_overlapped_bgcs = partially_overlapped_bgcs,
                                                               non_overlapped_bgcs = non_overlapped_bgcs)

        dereplicated_bgcs = output_shared_contigs['dereplicated_bgcs']
        overlapped_bgcs = output_shared_contigs['overlapped_bgcs']
        partially_overlapped_bgcs = output_shared_contigs['partially_overlapped_bgcs']
        non_overlapped_bgcs = output_shared_contigs['non_overlapped_bgcs']

        output_non_shared_contigs = dnsc.dereplicate_non_shared_contigs(metadata_df1= metadata_df1, 
                                                                        metadata_df2 = metadata_df2,
                                                                        dereplicated_bgcs = dereplicated_bgcs,
                                                                        non_overlapped_bgcs = non_overlapped_bgcs)

        dereplicated_bgcs = output_non_shared_contigs['dereplicated_bgcs']
        non_overlapped_bgcs = output_non_shared_contigs['non_overlapped_bgcs']
        
        ## Update metadata
        metadata = metadata[next_i: ]
        output_dir_dereplicated = f'{output_dir}/dereplicated'
        metadata_dereplicated_bgc = outputs.create_df(input_dict = dereplicated_bgcs, 
                                                      output_dir = output_dir_dereplicated,
                                                      df = True)
        metadata[0] = metadata_dereplicated_bgc
        
        if len(metadata) == 1:
            
            output_dict = {'dereplicated_bgcs': dereplicated_bgcs, 
                           'overlapped_bgcs': overlapped_bgcs, 
                           'partially_overlapped_bgcs': partially_overlapped_bgcs, 
                           'non_overlapped_bgcs': non_overlapped_bgcs}

        else:
            
                        
            output_dict = recursive_dereplication(metadata = metadata, 
                                    dereplicated_bgcs = dereplicated_bgcs,
                                    overlapped_bgcs = overlapped_bgcs,
                                    partially_overlapped_bgcs = partially_overlapped_bgcs,
                                    non_overlapped_bgcs = non_overlapped_bgcs)
                 
        return(output_dict)

###############################################################################
## 7. Execute dereplication
###############################################################################

output = recursive_dereplication(metadata = metadata_list, 
                                dereplicated_bgcs = dereplicated_bgcs,
                                overlapped_bgcs = overlapped_bgcs,
                                partially_overlapped_bgcs = partially_overlapped_bgcs,
                                non_overlapped_bgcs = non_overlapped_bgcs)


dereplicated_bgcs = output['dereplicated_bgcs']
overlapped_bgcs = output['overlapped_bgcs']
partially_overlapped_bgcs = output['partially_overlapped_bgcs']
non_overlapped_bgcs = output['non_overlapped_bgcs']

###############################################################################
## 8. Create output dirs
###############################################################################

output_dir_overlapped = f'{output_dir}/overlapped'
output_dir_partially_overlapped = f'{output_dir}/partially_overlapped'
output_dir_non_overlapped = f'{output_dir}/non_overlapped'
output_dir_dereplicated = f'{output_dir}/dereplicated'

if not os.path.exists(output_dir_overlapped):
    os.makedirs(output_dir_overlapped)

if not os.path.exists(output_dir_partially_overlapped):
    os.makedirs(output_dir_partially_overlapped)

if not os.path.exists(output_dir_non_overlapped):
    os.makedirs(output_dir_non_overlapped)

if not os.path.exists(output_dir_dereplicated):
    os.makedirs(output_dir_dereplicated)

###############################################################################
## 9. Sort data: create links
###############################################################################

# overlapped
outputs.create_links(input_dict = overlapped_bgcs, 
                     output_dir = output_dir_overlapped)

# partially overlapped
outputs.create_links(input_dict =  partially_overlapped_bgcs, 
                     output_dir = output_dir_partially_overlapped)

# non overlapped
outputs.create_links(input_dict = non_overlapped_bgcs, 
                     output_dir = output_dir_non_overlapped)

# dereplicated
outputs.create_links(input_dict = dereplicated_bgcs, 
                     output_dir = output_dir_dereplicated)

###############################################################################
## 10. Create DFs
###############################################################################

outputs.create_df(input_dict = dereplicated_bgcs,
                  output_dir = output_dir_dereplicated)

outputs.create_df(input_dict = overlapped_bgcs,
                  output_dir = output_dir_overlapped)

outputs.create_df(input_dict = non_overlapped_bgcs,
                  output_dir = output_dir_non_overlapped)

outputs.create_df(input_dict = partially_overlapped_bgcs,
                  output_dir = output_dir_partially_overlapped)
