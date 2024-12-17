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
sys.path.append(f'{current_dir}/src')
from bgc_derep import dsc, dnsc, outputs 

################################################################################
# 2. Parse parameters
################################################################################

parser = argparse.ArgumentParser(prog='dereplicate_bgcs.py', 
                                 description='Dereplicate BGCs annotated with different tools.')
parser.add_argument("--input_dir", default = None, help = "The input directory having the metadata tables (tsv) generated with the run_* modules.")
parser.add_argument("--overlap_thres", default = "0.75", 
                    help="Percentage of overlap (in relation to longest BGCs) to determine if two BGCs are overlapped or partially overlapped.")
parser.add_argument("--metadata", default = None, help="A comma-separated list of TSV tables containing the metadata of annotated BGCs. This parameter is ignored if --input_dir is used.")
parser.add_argument("--overwrite", action="store_true", help="Overwrite output directory.")
parser.add_argument("--output_dir", help = "The output directory where non-overlapped, partially overlapped and overlapped BGCs sequences and metadata tables (tsv) are saved.")

args = parser.parse_args()
input_dir = args.input_dir
overlap_thres = args.overlap_thres
metadata = args.metadata
overwrite = args.overwrite
output_dir = args.output_dir

################################################################################
# 3. Check input data is provided
################################################################################

if metadata is None and os.path.exists(input_dir) is False:
    print('Please provide metadata tables of annotated BGCs using either the --metadata or --input_dir arguments.')
    sys.exit(1)

################################################################################
# 4. Create or remove output directory
################################################################################

if output_dir is None:
    output_dir = os.path.join(os.path.dirname(os.path.abspath(input_dir)), "sorted")

if os.path.exists(output_dir) is True:
    
    if args.overwrite is False:
        print(f'Output dir {output_dir} already exists. Use --overwrite to overwrite')
        sys.exit(0)
    
    if args.overwrite is True:
        try:
            shutil.rmtree(output_dir)
        except OSError as e:
            print(f'Error: {e}')
        try:       
            os.makedirs(output_dir)
        except Exception as e:
            print(f'Error: {e}')
            sys.exit(1)
        
if os.path.exists(output_dir) is False:
    try:       
        os.makedirs(output_dir)
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)
        
################################################################################
# 4. Load the metadata tables and check if exist
################################################################################

if input_dir is not None:
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
        
        # dereplicate shared contigs
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
        
        # dereplicate non-shared contigs
        output_non_shared_contigs = dnsc.dereplicate_non_shared_contigs(metadata_df1= metadata_df1, 
                                                                        metadata_df2 = metadata_df2,
                                                                        dereplicated_bgcs = dereplicated_bgcs,
                                                                        non_overlapped_bgcs = non_overlapped_bgcs)

        dereplicated_bgcs = output_non_shared_contigs['dereplicated_bgcs']
        non_overlapped_bgcs = output_non_shared_contigs['non_overlapped_bgcs']
        
        ## Update metadata list
        metadata = metadata[next_i: ]
        output_dir_dereplicated = f'{output_dir}/dereplicated'
        output_dir_dereplicated_gbks = f'{output_dir}/dereplicated/gbks'
        metadata_dereplicated_bgc = outputs.create_df(input_dict = dereplicated_bgcs, 
                                                      output_dir = output_dir_dereplicated,
                                                      output_dir_gbks = output_dir_dereplicated_gbks,
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

# overlapped dir
output_dir_overlapped = f'{output_dir}/overlapped'
output_dir_overlapped_gbks = f'{output_dir}/overlapped/gbks'

if not os.path.exists(output_dir_overlapped):
    try:       
        os.makedirs(output_dir_overlapped_gbks)
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)

# partially overlapped dir
output_dir_partially_overlapped = f'{output_dir}/partially_overlapped'
output_dir_partially_overlapped_gbks = f'{output_dir}/partially_overlapped/gbks'

if not os.path.exists(output_dir_partially_overlapped):
    try:       
        os.makedirs(output_dir_partially_overlapped_gbks)
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)

# non overlapped dir   
output_dir_non_overlapped = f'{output_dir}/non_overlapped'
output_dir_non_overlapped_gbks = f'{output_dir}/non_overlapped/gbks'

if not os.path.exists(output_dir_non_overlapped):
    try:       
        os.makedirs(output_dir_non_overlapped_gbks)
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)

# dereplicated dir
output_dir_dereplicated = f'{output_dir}/dereplicated'
output_dir_dereplicated_gbks = f'{output_dir}/dereplicated/gbks'

if not os.path.exists(output_dir_dereplicated):
    try:       
        os.makedirs(output_dir_dereplicated_gbks)
    except Exception as e:
        print(f'Error: {e}')
        sys.exit(1)

###############################################################################
## 9. Sort data: create links
###############################################################################

# overlapped
outputs.create_links(input_dict = overlapped_bgcs, 
                     output_dir = output_dir_overlapped_gbks)

# partially overlapped
outputs.create_links(input_dict =  partially_overlapped_bgcs, 
                     output_dir = output_dir_partially_overlapped_gbks)

# non overlapped
outputs.create_links(input_dict = non_overlapped_bgcs, 
                     output_dir = output_dir_non_overlapped_gbks)

# dereplicated
outputs.create_links(input_dict = dereplicated_bgcs, 
                     output_dir = output_dir_dereplicated_gbks)

###############################################################################
## 10. Create DFs
###############################################################################

outputs.create_df(input_dict = dereplicated_bgcs,
                  output_dir = output_dir_dereplicated,
                  output_dir_gbks = output_dir_dereplicated_gbks)

outputs.create_df(input_dict = overlapped_bgcs,
                  output_dir = output_dir_overlapped,
                  output_dir_gbks = output_dir_overlapped_gbks)

outputs.create_df(input_dict = non_overlapped_bgcs,
                  output_dir = output_dir_non_overlapped,
                  output_dir_gbks = output_dir_non_overlapped_gbks)

outputs.create_df(input_dict = partially_overlapped_bgcs,
                  output_dir = output_dir_partially_overlapped,
                  output_dir_gbks = output_dir_partially_overlapped_gbks)

print("dereplicate executed successfully")