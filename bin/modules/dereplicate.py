#! /usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import pandas as pd
import os
import csv
import time
import argparse
import sys

current_dir = os.getcwd()
sys.path.append(f'{current_dir}/src')
from dereplicate import dnsc, dsc, outputs 

start_time = time.time()

################################################################################
# 2. Parse parameters
################################################################################

parser = argparse.ArgumentParser(prog='dereplicate_bgcs.py', description='Dereplicate BGCs annotated with different tools')
parser.add_argument("--overlap_thres", default = "0.75", 
                    help="Percentage of overlap (in relation to longest BGCs) to determine if two BGCs are overlapped or partially overlapped")
parser.add_argument("--metadata", help="Comma separated list of tsv tables containing the metadata of annotated BGCs")
parser.add_argument("--output_dir", help = "The output directory where non-overlapped, partially overlapped and overlapped BGCs sequences and metadata tables (tsv)")

args = parser.parse_args()
overlap_thres = args.overlap_thres
metadata = args.metadata
output_dir = args.output_dir

################################################################################
# 3. Load the metadata tables and check if exist
################################################################################

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
# 4. Initialize dictionaries
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
## 5. Define recursive function
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
        
        ## Updata metadata
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
## 6. Execute dereplication
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
## 7. Create output dirs
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
## 7. Sort data: create links
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
## 8. Create DFs
###############################################################################

outputs.create_df(input_dict = dereplicated_bgcs,
                  output_dir = output_dir_dereplicated)

outputs.create_df(input_dict = overlapped_bgcs,
                  output_dir = output_dir_overlapped)

outputs.create_df(input_dict = non_overlapped_bgcs,
                  output_dir = output_dir_non_overlapped)

outputs.create_df(input_dict = partially_overlapped_bgcs,
                  output_dir = output_dir_partially_overlapped)
