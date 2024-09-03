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
parser.add_argument("--input_bam", help="Input bam file.")
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

# BiG-SLICE parameters
parser.add_argument("--threshold_pct", default = 0.1, help="BiG-SLICE - Calculate clustering threshold (T) based on a random sampling of pairwise distances between the data, taking the N-th percentile value as the threshold.")

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

# MMSeqs taxonomy parameters
tax_lineage = args.tax_lineage 

# BiG-SLICE parameters
threshold_pct = args.threshold_pct

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

###############################################################################
## 6. Import precomputed BGC annotations
###############################################################################

# antismash
if antismash_output is not None and os.path.exists(antismash_output):
    antismash_output_current_dir = f"{bgc_annot_inter_dir}/antismash"
    antismash_output = os.path.abspath(antismash_output)
    try:
        shutil.copytree(antismash_output, antismash_output_current_dir)
    except OSError as e:
        print(f'Error: {e}')
    # set annotation flag: antismash was copied 
    antismash_annot_flag = True

# deepbgc
if deepbgc_output is not None and os.path.exists(deepbgc_output):
    deepbgc_output_current_dir = f"{bgc_annot_inter_dir}/deepbgc"
    deepbgc_output = os.path.abspath(deepbgc_output)
    try:
        shutil.copytree(deepbgc_output, deepbgc_output_current_dir)
    except OSError as e:
        print(f'Error: {e}')
    # set annotation flag: deepbgc was copied 
    deepbgc_annot_flag = True

# gecco 
if gecco_output is not None and os.path.exists(gecco_output):
    gecco_output_current_dir = f"{bgc_annot_inter_dir}/gecco"
    gecco_output = os.path.abspath(gecco_output)
    try:
        shutil.copytree(gecco_output, gecco_output_current_dir)
    except OSError as e:
        print(f'Error: {e}')
    # set annotation flag: gecco was copied 
    gecco_annot_flag = True

###############################################################################
## 7. BGC metadata
###############################################################################

###############################################################################
#### 7.1 anitsmash BGC metadata
###############################################################################

if antismash_annot_flag == True:

    # find all antismash GBK output files 
    antismash_output_gbks = find_files(input_dir = antismash_output_current_dir,
                                       pattern = ".*.region.*.gbk")
    
    try:
        os.makedirs(f'{antismash_output_current_dir}/gbks')
    except FileExistsError:
        print(f"Directory '{antismash_output_current_dir}/gbks' already exists.")
        sys.exit()

    for gbk in antismash_output_gbks:
            shutil.move(gbk, f'{antismash_output_current_dir}/gbks')

    antismash_annot_metadata = antismash_annot_parser(input_dir = f'{antismash_output_current_dir}/gbks',
                                                      sample_name = sample_name,
                                                      input_fasta = input_sample)

    if antismash_annot_metadata is not None:
        antismash_annot_metadata_tsv = f'{output_tables_dir}/antismash_annot_metadata.tsv'
        antismash_annot_metadata.to_csv(antismash_annot_metadata_tsv, sep='\t', index=False)
    else:
        antismash_annot_metadata_tsv = None

###############################################################################
#### 7.2 deepBGC BGC metadata
###############################################################################

if deepbgc_annot_flag == True:

    # find the sinlge GBK output file from deepBGC
    deepbgc_output_gbk = find_files(input_dir = deepbgc_output_current_dir, 
                                    pattern = ".bgc.gbk")
    if len(deepbgc_output_gbk) > 1:
        print("Error message:\nMore than one deepBGC GBK file output found\nThere should be only one.")
        print(deepbgc_output_gbk)
        sys.exit() 

    if len(deepbgc_output_gbk) == 1:
        # split multiple GBK into different files
        gbk_splitter(input_file = deepbgc_output_gbk[0],
                     output_dir = f'{deepbgc_output_current_dir}/gbks', 
                     sample_name = sample_name)

        # Filter deebBGC outputs with low score
        gbk_filter(input_dir = f'{deepbgc_output_current_dir}/gbks',
                   sample_name = sample_name,
                   deepbgc_score_thres = deepbgc_score_thres,
                   deepbgc_cds_count_thres = deepbgc_cds_count_thres,
                   output_dir = f'{deepbgc_output_current_dir}/gbks_filtered')

    deepbgc_annot_metadata = deepbgc_annot_parser(input_dir = f'{deepbgc_output_current_dir}/gbks',
                                                  sample_name = sample_name,
                                                  input_fasta = input_sample)


    if deepbgc_annot_metadata is not None:
        deepbgc_annot_metadata_tsv = f'{output_tables_dir}/deepbgc_annot_metadata.tsv'
        deepbgc_annot_metadata.to_csv(deepbgc_annot_metadata_tsv, sep='\t', index=False)
    else:
        deepbgc_annot_metadata_tsv = None

###############################################################################
#### 7.3 gecco BGC metadata
###############################################################################

if gecco_annot_flag == True:

    # copy all gecco output gbk files into a gbk folder and rename these following \
    # the format region*.gbk (as antiSMASH).
    gecco_gbk_list = find_files(input_dir =  gecco_output_current_dir, 
                                pattern = r'.cluster_[0-9]+.gbk')

    if not os.path.exists(f'{gecco_output_current_dir}/gbks'):
        os.makedirs(f'{gecco_output_current_dir}/gbks')

    for file_path in gecco_gbk_list:
        file = os.path.basename(file_path)
        pattern = r'_cluster_(\d+)\.gbk'
        match = re.search(pattern, file)
        number = int(match.group(1))
        file_renamed = re.sub(r'_cluster_\d+.gbk','.region{:03d}.gbk'.format(number), file)
        file_path_renamed = os.path.join(f'{gecco_output_current_dir}/gbks', file_renamed)
        shutil.copy(file_path, file_path_renamed)

    gecco_annot_metadata = gecco_annot_parser(input_dir = f'{gecco_output_current_dir}',
                                              sample_name = sample_name,
                                              input_fasta = input_sample)

    if gecco_annot_metadata is not None:
        gecco_annot_metadata_tsv = f'{output_tables_dir}/gecco_annot_metadata.tsv'
        gecco_annot_metadata.to_csv(gecco_annot_metadata_tsv, sep='\t', index=False)
    else:
        gecco_annot_metadata_tsv = None

###############################################################################
#### 8. Format GBKs
###############################################################################

###############################################################################
#### 8.1. Format deepBGC GBKs
###############################################################################

if deepbgc_annot_metadata_tsv is not None:

    format_gbks(input_dir = f'{deepbgc_output_current_dir}/gbks',
                input_tsv = deepbgc_annot_metadata_tsv,
                output_dir = f'{deepbgc_output_current_dir}/gbks',
                tool = "deepbgc")

    format_gbks(input_dir = f'{deepbgc_output_current_dir}/gbks_filtered',
                input_tsv = deepbgc_annot_metadata_tsv,
                output_dir = f'{deepbgc_output_current_dir}/gbks_filtered',
                tool = "deepbgc")

###############################################################################
#### 8.2. Format gecco GBKs 
###############################################################################

if gecco_annot_metadata_tsv is not None:

    format_gbks(input_dir = f'{gecco_output_current_dir}',
                input_tsv = gecco_annot_metadata_tsv,
                output_dir = f'{gecco_output_current_dir}/gbks',
                tool = "gecco")

###############################################################################
#### 9. Dereplicate
###############################################################################

###############################################################################
#### 8.1 Dereplicate antisSMASH vs deepbgc
###############################################################################

dereplicated_annot_metadata_1_tsv =  f'{output_tables_dir}/dereplicated_annot_metadata_1.tsv'

dereplicate_bgc(metadata_tsv1 = antismash_annot_metadata_tsv,
                metadata_tsv2 = deepbgc_annot_metadata_tsv,
                sample_name = sample_name,
                output_dir = f'{bgc_annot_sorted_dir}/tmp1',
                output_tsv = dereplicated_annot_metadata_1_tsv)

###############################################################################
#### 8.2 Dereplicate: dereplicated vs gecco
###############################################################################

if not os.path.exists(dereplicated_annot_metadata_1_tsv):
    dereplicated_annot_metadata_1_tsv = None

dereplicated_annot_metadata_2_tsv =  f'{output_tables_dir}/dereplicated_annot_metadata_2.tsv'

dereplicate_bgc(metadata_tsv1 = dereplicated_annot_metadata_1_tsv,
                metadata_tsv2 = gecco_annot_metadata_tsv,
                sample_name = sample_name,
                output_dir = f'{bgc_annot_sorted_dir}',
                output_tsv = dereplicated_annot_metadata_2_tsv)

###############################################################################
#### 8.3 Dereplicate: Remove tmp files
###############################################################################

if os.path.exists(f'{bgc_annot_sorted_dir}/tmp1'):
    shutil.rmtree(f'{bgc_annot_sorted_dir}/tmp1')

###############################################################################
#### 4.2. Estimate BGC class coverage
###############################################################################

# bgc_coverage_df = get_coverage(input_dir=f'{bgc_annot_sorted_dir}/dereplicated/{sample_name}', 
#                                input_bam=input_bam, 
#                                sample_name=sample_name)

# tsv_file_path = f"{output_tables_dir}/bgc_abund.tsv"
# bgc_coverage_df.to_csv(tsv_file_path, sep='\t', index=False)

# ###############################################################################
# ## 5. Perform taxonomic annotation 
# ###############################################################################

# ###############################################################################
# #### 5.1. Gather all contigs in which a BGC was annotated 
# ###############################################################################

# bgc_acc_list = bgc_coverage_df["acc"].tolist()
# bgc_contigs = extract_sequences_by_ids(fasta_file = input_sample, 
#                                        target_ids = bgc_acc_list)

# records = []
# for seq_id, sequence in bgc_contigs.items():
#   record = f">{seq_id}\n{sequence}\n"
#   records.append(record)

# bgc_contigs_file = f'{bgc_taxa_dir}/bgc_contigs.fasta'

# with open(bgc_contigs_file, 'w') as output_handle:
#   output_handle.write("".join(records))

# ###############################################################################
# #### 5.2. Run MMSeqs taxonomy
# ###############################################################################

# command_mmseqs = f"{current_directory}/src/run_mmseqs_taxonomy.sh \
#   {bgc_taxa_dir}/bgc_contigs.fasta \
#   {bgc_taxa_dir}/bgc_taxonomy \
#   --threads {threads} \
#   --tax-lineage {tax_lineage} \
#   -v 3"

# command_mmseqs = subprocess.run(command_mmseqs, 
#                                 shell=True, 
#                                 stdout=subprocess.PIPE, 
#                                 stderr=subprocess.PIPE, 
#                                 text=True)

# if command_mmseqs.returncode == 0:
#     print("MMSeqs taxonomy executed successfully")
# else:
#     print("Error executing MMSeqs taxonomy")
#     print("Error message:\n", command_mmseqs.stderr)

# ###############################################################################
# #### 5.3. Format MMSeqs taxonomy output
# ###############################################################################

# bgc_taxonomy_file = f"{bgc_taxa_dir}/bgc_taxonomy_lca.tsv"
# bgc_taxonomy_file_df = pd.read_csv(bgc_taxonomy_file, sep='\t', header=None)

# column_names = ['acc', 'taxid', 'rank', 'taxon_name', 'fragments_retained', \
#                 'fragments_annotated', 'fragments_consistent', 'support', 'lineage']
# bgc_taxonomy_file_df.columns = column_names

# bgc_taxonomy_file_formatted_df = bgc_taxonomy_file_df[["acc","taxon_name","lineage"]] 
# bgc_taxonomy_file_formatted_df.loc[:,["sample"]] = bgc_taxonomy_file_formatted_df["acc"].str.replace("__k[0-9].*", "", regex=True)
# bgc_taxonomy_file_formatted_df = bgc_taxonomy_file_formatted_df.loc[:,["sample", "acc", "taxon_name", "lineage"]]

# tsv_file_path = f"{output_tables_dir}/bgc_taxa.tsv"
# bgc_taxonomy_file_formatted_df.to_csv(tsv_file_path, sep='\t', index=False)

# ###############################################################################
# ## 6. BGC clustering
# ###############################################################################

# ###############################################################################
# #### 6.1 Create dataset.tsv file and taxonomy folder and files
# ###############################################################################

# # this path will change after doing the dereplication
# bigslice_input_dir = bgc_annot_inter_antismash_dir

# create_taxonomy_tables(bigslice_input_dir)

# create_dataset_table(bigslice_input_dir)

# ###############################################################################
# ### 6.2 Run BiG-SLICE querying mode
# ###############################################################################

# mibig_gcf_models_dir = f"{current_directory}/resources/mibig_gcf_models"

# command_bigslice = f"{current_directory}/src/run_bigslice.sh \
#    query {bigslice_input_dir} {mibig_gcf_models_dir} \
#   --num_threads {threads} \
#   --threshold_pct {threshold_pct} \
#   --query_name {sample_name}"

# command_bigslice = subprocess.run(command_bigslice, 
#                                 shell=True, 
#                                 stdout=subprocess.PIPE, 
#                                 stderr=subprocess.PIPE, 
#                                 text=True)

# if command_bigslice.returncode == 0:
#     print("BiG-SLICE executed successfully")
# else:
#     print("Error executing BiG-SLICE")
#     print("Error message:\n", command_bigslice.stderr)

# ###############################################################################
# ### 6.3. Format BiG-SLICE querying mode output
# ###############################################################################

# reports_db = f"{mibig_gcf_models_dir}/reports/1/data.db"
# result_db = f"{mibig_gcf_models_dir}/result/data.db"

# reports_con = sqlite3.connect(reports_db)
# reports_cur = reports_con.cursor()

# result_con = sqlite3.connect(result_db)
# result_cur = result_con.cursor()

# res = reports_cur.execute("SELECT gcf_membership.bgc_id, bgc_class.chem_subclass_id \
#                             FROM gcf_membership, bgc_class \
#                             WHERE gcf_membership.bgc_id=bgc_class.bgc_id")


# subclass_ids = [list(row)[1] for row in res]
# subclass_ids_uniq = list(dict.fromkeys(subclass_ids))

# query = "SELECT * FROM chem_subclass_map WHERE subclass_id IN ({})".format(', '.join(['?' for _ in subclass_ids_uniq]))

# result_cur.execute(query, subclass_ids_uniq)
# rows = result_cur.fetchall()

# print(rows)

# # print((res.fetchall()))

# reports_con.close()
# result_con.close()



