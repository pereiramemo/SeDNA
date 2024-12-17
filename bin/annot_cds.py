#! /usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import pandas as pd
import subprocess
import argparse
import os
import sys
import shutil
import pyhmmer

current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{current_dir}/src')
from bgc_annot import utilities
                                                       
################################################################################
# 2. Parse parameters
################################################################################

parser = argparse.ArgumentParser(prog='cds_annot.py', \
                                 description='Annotate CDS sequences with pyHMMER.')


parser.add_argument("--input_dir", help="Input directory with GenBank files.")
parser.add_argument("--pattern", default = ".*.gbk", help="Pattern in GenBank file names.")
parser.add_argument("--cds_id", default = "locus_tag", help="Key to be used as CDS ID.")
parser.add_argument("--create_fasta", default = True, action=argparse.BooleanOptionalAction, help="True or False to create fasta files from CDSs.")
parser.add_argument("--annotate", default = True, action=argparse.BooleanOptionalAction, help="True or False to perform the CDSs annotation using hmmsearch.")
parser.add_argument("--resolve_doms", default = True, action=argparse.BooleanOptionalAction, help="Resolve domain structures using cath-resolve-hits.")
parser.add_argument("--overwrite", default = False, action=argparse.BooleanOptionalAction, help="Overwrite output directory.")
parser.add_argument("--hmms_db", default = f'{resources_dir}/databases/Pfam-A.hmm', help="Path to HMMs database.")
parser.add_argument("--evalue_thres", default = 1e-3, help="hmmsearch - e-value threshold.")
parser.add_argument("--cut_ga", default = True, action=argparse.BooleanOptionalAction, help="hmmsearch - Use the gathering bitscores for sequence inclusion")
parser.add_argument("--threads", default = 4, help="hmmsearch - Number of threads")
parser.add_argument("--output_dir", help="Output directory.")

args = parser.parse_args()
input_dir = args.input_dir
pattern = args.pattern
cds_id = args.cds_id
create_fasta = args.create_fasta
annotate = args.annotate
resolve_doms = args.resolve_doms
overwrite = args.overwrite
hmms_db = args.hmms_db
evalue_thres = float(args.evalue_thres)
threads = int(args.threads)
cut_ga = args.cut_ga
output_dir = args.output_dir

################################################################################
# 3. Create output dir
################################################################################

if output_dir is None:
    output_dir = os.path.join(input_dir, "cds_annot")

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
# 4. Check input dir
################################################################################

if not os.path.exists(input_dir):
    print(f'Input dir {input_dir} does not exist')
    sys.exit(1)

################################################################################
# 5. Create fasta files
################################################################################

path2gbks = utilities.find_files(input_dir = input_dir, pattern = pattern)

if create_fasta == True:
    paths2gbk = utilities.find_files(input_dir = input_dir,  
                                     pattern = pattern)

    if len(paths2gbk) == 0:
        print("No gbk file found")
        sys.exit()

    path2names = []

    for gbk in paths2gbk:
        file_name = os.path.basename(os.path.splitext(gbk)[0])
        file_path = os.path.join(output_dir, file_name)
        fasta = os.path.join(output_dir, f'{file_name}.fasta')
        path2names.append(file_path)
        utilities.extract_cds(input_gbk = gbk, 
                              cds_id = cds_id, 
                              output_fasta = fasta)

################################################################################
# 6. HMM annotation
################################################################################

if annotate == True:

    for file in path2names:

        # IO files
        fasta = f'{file}.fasta' 
        domtblout = f'{file}_domtblout.tsv'
        hmmout = f'{file}_hmmout.tsv'  

        if cut_ga is not True:
            with pyhmmer.plan7.HMMFile(hmms_db) as hmms:
                with pyhmmer.easel.SequenceFile(fasta, digital=True) as seqs:
                    with open(domtblout, "wb") as output_file:
                        for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, domE = evalue_thres, cpus = threads):
                            hits.write(output_file, format = "domains", header = False)
        else:
            with pyhmmer.plan7.HMMFile(hmms_db) as hmms:
                with pyhmmer.easel.SequenceFile(fasta, digital=True) as seqs:
                    with open(domtblout, "wb") as output_file:
                        for hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, domE = evalue_thres, 
                                                            bit_cutoffs="gathering", cpus = threads):
                            hits.write(output_file, format = "domains", header = False)
                                            
################################################################################
# 7. Resolve domain structure
################################################################################

if resolve_doms == True:

    for file in path2names:

        domtblout = f'{file}_domtblout.tsv'
        output_tsv = f'{file}_annotdoms.tsv'
        output_resolved_tsv = f'{file}_annotdoms_resolved.tsv'

        command_cath_resolve_hits = f"cath-resolve-hits \
                                --input-format hmmer_domtblout \
                                {domtblout} > {output_resolved_tsv}"

        result_cath_resolve_hits = subprocess.run(command_cath_resolve_hits, 
                                          shell=True, 
                                          stdout=subprocess.PIPE, 
                                          stderr=subprocess.PIPE, 
                                          text=True)

        if result_cath_resolve_hits.returncode != 0:
            print("Error executing cath-resolve-hits")
            print("Error message:\n", result_cath_resolve_hits.stderr)
            sys.exit()

print("annot_cds executed successfully")