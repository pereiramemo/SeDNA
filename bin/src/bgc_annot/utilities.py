###############################################################################
# 1. find_files
###############################################################################

import os
import re

def find_files(input_dir : str = None, pattern : str = r".*\.region\d+.gbk") -> list:

    '''
    input_dir:
      type: str
      Directory where the search is performed.
    pattern:
        Type: str
        Pattern of files names to search.  
    returns:
        A list of paths to file names found.
   '''
    
    matching_files = []
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if re.search(pattern, file):
                matching_files.append(os.path.join(root, file))
    return matching_files

'''
Description
Given a directory and pattern, the function recursively finds all the files that match the such a pattern (i.e., coded as a regular expression), and returns a list with all the paths to file names found.
'''

###############################################################################
# 2. find_dirs
###############################################################################

import os
import re

def find_dirs(input_dir : str = None, pattern : str = r".*\.region\d+.gbk") -> list:

    '''
    input_dir:
      type: str
      Directory where the search is performed.
    pattern:
        Type: str
        Pattern of dirs names to search.  
    returns:
        A list of paths to dirs names found.
   '''

    matching_dirs = []
    for root, dirs, files in os.walk(input_dir):
        for dir in dirs:
            if re.search(pattern, dir):
                matching_dirs.append(os.path.join(root, dir))
    return matching_dirs
    
# find_dirs(input_dir ='/home/jovyan/sandbox/development/epereira/bioprospecting//results/oceandna_redu/bgc_annt/deepbgc', pattern = "OceanDNA-.*\d+")

'''
Description
Given a directory and pattern, the function recursively finds all the dirs that match the such a pattern (i.e., coded as a regular expression), and returns a list with all the paths to dir names found.
'''

###############################################################################
# 3. get_features
###############################################################################

from Bio import SeqIO

def get_features(record : str = None, feature_name : str = None, qualifier_name : str = None) -> str:
    '''
    record:
        A GBK object generated with Bio.
    feature_name:
        Name of the feature to parse.
    qualifier_name:
        Name of the qualifier to parse within feature_name.
    returns:
        The value of qualifier_name within feature name.
    '''
    for feature in record.features:
        if feature.type == feature_name:
            return feature.qualifiers[qualifier_name][0]      
        
'''
Description
This function takes a GBK object as input and extracts the qualifier `qualifier_name` within `feature_name`.
This is a simple function that is used in other functions defined below, such as `antismash_parse_gbk`, and `deepbgc_parse_gbk`.
'''

###############################################################################
# 4. get_feature_location
###############################################################################

from Bio import SeqIO   

def get_feature_location(record : str = None, feature_name : str = None) -> list:
    '''
    record:
        A GBK object generated with Bio.
    feature_name:
        Name of the feature to parse.
    returns:
        A list containing the start and end coordinates for feature_name.
    '''
    for feature in record.features:
        if feature.type == feature_name:
            return [int(feature.location.start), int(feature.location.end)]

'''
Description
This function extracts the start and end coordinates of `feature_name`, and returns these as a two value list (i.e., [start,end]).
This is a simple function that is used in other functions defined below, such as `antismash_parse_gbk`, and `get_coverage`.
'''

###############################################################################
# 5. get_seqs_length
###############################################################################

from Bio import SeqIO
import gzip
import os
            
def get_seqs_length(input_fasta: str = None) -> dict:

    '''
    input_fasta:
    Path to a fasta file.
    returns:
    A dictionary containing the sequence id as key and the sequence length as value.
    '''
    def is_gzipped(file_path):
        extension = os.path.splitext(file_path)[1]
        if extension == ".gz":
            return True
        else:
            return False
    sequences_lengths = {}

    if is_gzipped(input_fasta):
        with gzip.open(input_fasta, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences_lengths[record.id] = len(record.seq)
            return sequences_lengths
    else:
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequences_lengths[record.id] = len(record.seq)
        return sequences_lengths
    
'''
Description
This is a simple function to determine the length of each sequence in a fasta file. This is a simple function that is used in other functions defined below, such as `antismash_parse_gbk`, and `deepbgc_parse_gbk`.
'''

###############################################################################
# 6. gbk_splitter
###############################################################################

import sys
import os
import shutil
import re
from Bio import SeqIO 

def gbk_splitter(input_file : str = None, output_dir : str = None, 
                 sample_name : str = None): 
  
    '''
    input_file
        Single GBK file containing multiple entries.
    output_dir
        Output dir where all splitted GBK files are saved.
    '''

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    try:
        os.makedirs(output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
  
    # regions counter dict
    bgc_region_names = {}
    bgc_contig_id_counter = {}

    # open deepbgc output GBK file
    for gbk in SeqIO.parse(input_file, "gb"):
        # acc
        acc = gbk.annotations['accessions'][0]
        # contig_id
        pattern = re.compile('_[0-9]+-[0-9]+$')
        contig_id = pattern.sub('', acc)
        # create dict
        if contig_id not in bgc_contig_id_counter:
            bgc_contig_id_counter[contig_id] = 1
            bgc_region_names[acc] =  "{0}.region{1:03d}.gbk".format(contig_id,
                                                                  bgc_contig_id_counter[contig_id])
        else:
            bgc_contig_id_counter[contig_id] += 1
            bgc_region_names[acc] =  "{0}.region{1:03d}.gbk".format(contig_id,
                                                                  bgc_contig_id_counter[contig_id])


        output_file = os.path.join(output_dir, bgc_region_names[acc])
        SeqIO.write(gbk, output_file, "gb")
                
'''
Description
This function takes a single GBK file containing multiple entries (concatenated) and splits it into multiple GBK, one per sequence, which are saved in the `output_dir/gbks` folder.
The names of the GBK files follow the nomenclature from antiSMASH, that is, the contig_id followed by a counter. 
The purpose of this function is to have single entry GBKs to be able to lead these with the SeqIO module.
'''

###############################################################################
# 7. gbk_filter
###############################################################################

import sys
import os
import shutil
import re
from Bio import SeqIO 

def gbk_filter(input_dir : str = None, sample_name : str = None,
               deepbgc_score_thres : float = 0.75, deepbgc_cds_count_thres: int = 2, 
               output_dir : str = None): 
  
    '''
    input_dir:
        Directory containing splitted GBK files to be removed.
    score_thres:
        DeepBGC score threshold to filter GBK files.
    cds_count:
        CDS count threshold to filter GBK files.
    sample_name:
        Name of the sample.
    '''
    
    '''
    Start subfinctions def
    '''
    def count_cds(record):
        cds_count = 0
        for feature in record.features:
            if feature.type == "CDS":
                cds_count += 1
        return int(cds_count)
    '''
    End subfinctions def
    '''

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    try:
        os.makedirs(output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

    paths2gbk = find_files(input_dir = input_dir, 
                        pattern = f"{sample_name}.*region.*.gbk")
    
    for input_gbk in paths2gbk:
        with open(input_gbk, "r") as gbk_handle:
            record = SeqIO.read(gbk_handle, "genbank")
            
            # get deepbgc score
            deepbgc_score = get_features(record = record, 
                                         feature_name = "cluster",
                                         qualifier_name = "deepbgc_score")
            deepbgc_score = float(deepbgc_score)
            
            # count number of CDS
            cds_count = count_cds(record = record)
            
            # filter 
            if deepbgc_score < deepbgc_score_thres or cds_count < deepbgc_cds_count_thres:
                shutil.move(input_gbk, output_dir)

'''
Description

'''

###############################################################################
# 8. antismash_annot_parser
###############################################################################

import pandas as pd
from Bio import SeqIO

def antismash_annot_parser(input_dir: str = None, sample_name: str = None, 
                           tool : str = "antismash", input_fasta: str = None) -> pd.DataFrame:
  
    '''
    input_dir:
        The input folder containing antismash subfolders.
    input_fasta:
        The input fasta file (assembly) in which the BGC sequences were annotated     
    returns: 
        A data frame containing the sample name, BGC id (i.e., acc), start and end coordinates.
    '''

    '''
    Start subfuncions definitions
    '''
    
    def gbk_parser(paths2gbk : list = None, fasta_lengths_dict : dict = None):  
        acc2x = dict()   
        i = 1
        for input_gbk in paths2gbk:
            with open(input_gbk, "r") as gbk_handle:
                record = SeqIO.read(gbk_handle, "genbank")
                # acc
                acc = record.annotations['accessions'][0]
                # contig_id
                contig_id = record.name
                # contig_length
                contig_length = fasta_lengths_dict[contig_id]
                # bgc class
                bgc_class = get_features(record,"cand_cluster", "product")
                # start and end coordinates and length
                location = get_feature_location(record,"cand_cluster")
                start = int(location[0])
                end = int(location[1])
                length = end - start # 1 is not added given that coords are zero based
                # on edge
                on_edge = get_features(record,"cand_cluster", "contig_edge")  
                # input gbk
                input_gbk = os.path.abspath(input_gbk)
                # add index
                acc_i = acc + "-" + str(i)
                i += 1
                # create output dict
                acc2x[acc_i] = {"acc": acc, "bgc_class":bgc_class, \
                                "start":start, "end":end, "length":length, \
                                "on_edge": on_edge, "contig_id": contig_id, \
                                "contig_length": contig_length, 'file':input_gbk, \
                                "sample_name": sample_name, 'tool': tool}

        return(acc2x)

    '''
    End subfuncions definitions
    '''
    # find all antismash BGCs
    antismash_gbks = find_files(input_dir = input_dir, pattern = r'.region[0-9]+.gbk')
    if len(antismash_gbks) == 0:
        return

    # get sequences length
    fasta_lengths_dict = get_seqs_length(input_fasta)
  
    # get all antismash BGC ids and coords
    antismash_gbks_dict = gbk_parser(antismash_gbks, fasta_lengths_dict)
   
    # convert to df
    antismash_gbks_df = pd.DataFrame.from_dict(antismash_gbks_dict, orient='index')
      
    return(antismash_gbks_df)

'''
Description
The function detects all antiSMASH GBK files (region[0-9]+.gbk) an parses each file to create a DF with the following columns: `acc`, `bgc_class`, `start`, `end`, `length`, `on_edge`, `contig_id`, `contig_length`, `file`, `sample_name`, and `tool` for each identified BGC sequence.
In addition, the assembled sequences, in which the BGC sequences were annotated, must be provided as an input file in order to determine the contig length. 
Note that since the acc provided by antiSMASH for each BGC is the same as the contig name, the unique identifiers of the BGCs are the rownames of the DF (acc + a counter). Otherwise, the function would only output the last BGCs of all BGCs found in the same contig.
'''

###############################################################################
# 9. deepbgc_annot_parser
###############################################################################

from Bio import SeqIO
import pandas as pd
import re

def deepbgc_annot_parser(input_dir: str = None, sample_name: str = None, 
                         tool : str = "deepbgc", input_fasta: str = None) -> pd.DataFrame:
    '''
    input_dir:
      The input folder containing deepbgc subfolders.
    input_fasta:
      The input fasta file (assembly) in which the BGC sequences were annotated.     
    returns: 
      A data frame containing the sample name, BGC id (i.e., acc), start and end coordinates.
    '''
    
    '''
    Start subfuncions definitions
    '''              
    def gbk_parser(paths2gbk : list = None, fasta_lengths_dict : dict = None):
        acc2x = dict()   
        i = 1
        for input_gbk in paths2gbk:
            with open(input_gbk, "r") as gbk_handle:
                record = SeqIO.read(gbk_handle, "genbank")
                # acc
                acc = record.annotations['accessions'][0]
                # contig_id
                pattern = re.compile('_[0-9]+-[0-9]+$')
                contig_id = pattern.sub('', acc)
                # start and end coords, and length
                coords = re.split(r'[-,_]',acc)
                start = int(coords[-2])
                end = int(coords[-1])
                length = end - start # 1 is not added given that coords are zero based
                # add index
                acc_i = acc + "-" + str(i)
                i += 1
                # bgc class
                bgc_class = get_features(record, "cluster", "product_class_score")
                bgc_class = bgc_class.replace(" ", "")
                bgc_class_list = bgc_class.split(',')
                bgc_class_dict = {key_value.split('=')[0]: float(key_value.split('=')[1]) for key_value in bgc_class_list}
                bgc_class = max(bgc_class_dict, key=lambda k: bgc_class_dict[k])
                # contig_length
                contig_length = fasta_lengths_dict[contig_id]
                # on edge
                if start < 3 or (contig_length - end) <= 30:
                    on_edge = "True"
                else:
                    on_edge = "False"
                # input gbk
                input_gbk = os.path.abspath(input_gbk)
                # create output dict
                acc2x[acc_i] = {"acc": acc, "bgc_class":bgc_class, \
                              "start":start, "end":end, "length":length, \
                              "on_edge": on_edge, "contig_id": contig_id, \
                              'contig_length': contig_length, 'file':input_gbk, \
                              "sample_name": sample_name, 'tool': tool}
        return(acc2x)

    '''
    End subfuncions definitions
    '''
 
    # find all deepbgc BGCs
    deepbgc_gbks = find_files(input_dir = input_dir, pattern = r'.region[0-9]+.gbk')
    if len(deepbgc_gbks) == 0:
        return
    
    # get sequences length 
    fasta_lengths_dict = get_seqs_length(input_fasta)
    
    # get all deepbgc BGC ids and coords
    deepbgc_gbks_dict = gbk_parser(deepbgc_gbks, fasta_lengths_dict)
   
    # convert to df
    deepbgc_gbks_df = pd.DataFrame.from_dict(deepbgc_gbks_dict, orient='index')
      
    return(deepbgc_gbks_df)

'''
Description
This function was designed to produce the same output as  `antismash_annot_parser()`. That is, creating a DF with the following columns: `acc`, `bgc_class`, `start`, `end`, `length`, `on_edge`, `contig_id`, `contig_length`, `file`, `sample_name`, and `tool` for each identified BGC sequence.
Since the input GBK files are generated by deepBGC, the field `on_edge` has to be determined based on the BGC coordinates and the contig length.  Thus, the assembled sequences, in which the BGC sequences were annotated, must be provided as an input file. The criteria to determine the on_edge field is: if the start position is less than 3 or if the end coordinate is less than 30 bp from the end of the contig. To determine the BGC class, the `produc_class` with the highest score is selected.
This function was designed to produce the same output as  `antismash_gbk_parser()`. That is, creating a DF with the following columns: `acc`, `bgc_class`, `start`, `end`, `length`, `on_edge`, `contig_id`, `contig_length`, `file`, and `sample_name`, for each identified BGC sequence.
Since the input GBK files are generated by deepBGC, and the fields `on_edge` has to be determined based on the BGC coordinates and the contig length.  Thus, the assembled sequences, in which the BGC sequences were annotated, must be provided as an input file. The criteria to determine the on_edge field is: if the start position is less than 3 or if the end coordinate is less than 30 bp from the end of the contig. To determine the BGC class, the `produc_class` with the highest score is selected.
'''

###############################################################################
# 10. gecco_annot_parser
###############################################################################

import pandas as pd
import re

def gecco_annot_parser(input_dir: str = None, sample_name: str = None, 
                       tool : str = "gecco", input_fasta: str = None) -> pd.DataFrame:
    '''
    input_dir:
      The input folder containing gecco subfolders.
    input_fasta:
      The input fasta file (assembly) in which the BGC sequences were annotated.     
    returns: 
      A data frame containing the sample name, BGC id (i.e., acc), start and end coordinates.
    '''
    
    '''
    Start subfuncions definitions
    '''              
    def tbl_parser(path2tsv : str = None, fasta_lengths_dict : dict = None):
        acc2x = dict()   
        i = 1
        df = pd.read_csv(path2tsv, sep='\t')  
        for index, row in df.iterrows():
            # acc
            acc = row['cluster_id']
            # contig_id
            contig_id = row['sequence_id']
            # start and end coords, and length
            start = int(row['start'])
            end = int(row['end'])
            length = end - start # 1 is not added given that coords are zero based        
            # add index
            acc_i = acc + "-" + str(i)
            i += 1
            # bgc class
            bgc_class = row['type']
            # contig_length
            contig_length = fasta_lengths_dict[contig_id]
            # on edge
            if start < 3 or (contig_length - end) <= 30:
                on_edge = "True"
            else:
                on_edge = "False"
            # input gbk
            pattern = r'_cluster_(\d+)'
            match = re.search(pattern, acc)
            extracted_number = int(match.group(1))
            acc_renamed = re.sub(r'_cluster_\d+','.region{:03d}'.format(extracted_number), acc)
            input_gbk_list = find_files(input_dir = input_dir, pattern = f'{acc_renamed}.gbk')
            input_gbk = input_gbk_list[0]
            # create output dict
            acc2x[acc_i] = {"acc": acc, "bgc_class":bgc_class, \
                            "start":start, "end":end, "length":length, \
                            "on_edge": on_edge, "contig_id": contig_id, \
                            'contig_length': contig_length, 'file':input_gbk, \
                            "sample_name": sample_name, 'tool': tool}
        return(acc2x)

    '''
    End subfuncions definitions
    '''
 
    # find all deepbgc BGCs
    gecco_tsv = find_files(input_dir = input_dir, pattern = '.clusters.tsv')
    if len(gecco_tsv) == 0:
        return
    
    # get sequences length 
    fasta_lengths_dict = get_seqs_length(input_fasta)
    
    # get all gecco BGC ids and coords
    gecco_tbl_dict = tbl_parser(gecco_tsv[0], fasta_lengths_dict)
   
    # convert to df
    gecco_tbl_df = pd.DataFrame.from_dict(gecco_tbl_dict, orient='index')
      
    return(gecco_tbl_df)

'''
Description
This function was designed to produce the same output as  `antismash_annot_parser()` and `deepbgc_annot_parser()`. That is, creating a DF with the following columns: `acc`, `bgc_class`, `start`, `end`, `length`, `on_edge`, `contig_id`, `contig_length`, `file`, `sample_name`, and `tool` for each identified BGC sequence.
Since `gecco` already generates a table (i.e., `*.clusters.tsv`) where several of these fields are included (also the `bgc_class`), it is not necessary to parse GBK files. However, similarly as performed in `deepbgc_annot_pareser()`, the field `on_edge` has to be determined based on the BGC coordinates and the contig length. Thus, the assembled sequences, in which the BGC sequences were annotated, must be provided as an input file. The criteria to determine the `on_edge` field is: if the start position is less than 3 or if the end coordinate is less than 30 bp from the end of the contig. 
A particular behavior of the function is that the `file` field is added (with the `find_files` function) utilizing files renamed following antiSMASH naming (i.e., *.region*.gbk), and copied into a `gbks` folder inside `gecco`’s output. Thus, this folder, with GBK files following such naming convention, has to be generated before running this function.
'''

###############################################################################
# 11. format_gbks
###############################################################################

# Adapted from Satria A. Kautsar
# Wageningen University & Research
# Bioinformatics Group
# Copyright (C) 2020

"""generate custom antiSMASH regiongbks"""

import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import pandas as pd

def format_gbks(input_dir: str = None, input_tsv: str = None,
                tool : str = None, output_dir: str = None):
  
    '''
    input_dir:
      Directory containing GBK files.
    input_csv: 
      A tsv file with three columns: the BGC id, and the start and end coordinates. 
    output_dir:
      The output directory.  
    returns:
      An output directory where all formatted GBKs, one per BGC sequences, following the naming and standards of antiSMASH.
    '''

    ############################################################################
    # Create output dir
    ############################################################################
  
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
  
    ############################################################################
    # Parse coordinates file: start and end positions
    ############################################################################
  
    bgc_metadata = {}
    df = pd.read_csv(input_tsv, sep='\t')
    for idx, line in df.iterrows():
    
        acc = line['acc']
        start = line['start']
        end = line['end']
        bgc_class = line['bgc_class']
        on_edge = line['on_edge']
        start = int(start)
        end = int(end)
    
        if acc not in bgc_metadata:
            bgc_metadata[acc] = []

        bgc_metadata[acc].append({
        "start": start,
        "end": end,
        "bgc_class": bgc_class,
        "on_edge": on_edge,
        "num": idx
        })

    ############################################################################
    # Create GBK 
    ############################################################################
  
    # find all deepbgc BGCs
    paths2gbk = find_files(input_dir = input_dir, pattern = r'.region[0-9]+.gbk')
  
    for input_gbk in paths2gbk:
        with open(input_gbk, "r") as gbk_handle:
      
            gbk = SeqIO.read(gbk_handle, "genbank")
            acc = gbk.annotations['accessions'][0]
      
            if tool == "deepbgc":
                contig_id = re.sub(r'_[0-9]+-[0-9]+$','',acc)
                gbk.name = contig_id

            if tool == "gecco":
                contig_id = re.sub(r'_cluster_[0-9]+$','',acc)
                gbk.name = contig_id
    
            if acc in bgc_metadata:
                # generate regiongbks
                for region_meta in bgc_metadata[acc]:
                    # add molecule type
                    gbk.annotations["molecule_type"] = "dna"
                    # make antiSMASH5-like structured comment
                    gbk.annotations["structured_comment"] = {
                    "antiSMASH-Data": {
                    "Version": "Simulated",
                    "Orig. start": region_meta["start"],
                    "Orig. end": region_meta["end"],
                    "Note": ("This is a simulated antiSMASH regiongbk!!!")}}

                    new_features = []

                    # make protocluster feature (only for detecting type)
                    proto_feature = SeqFeature(FeatureLocation(0, len(gbk)), type="protocluster")
                    new_features.append(proto_feature)
                    # make region feature
                    region_feature = SeqFeature(FeatureLocation(0, len(gbk)), type="region")
                    region_feature.qualifiers["contig_edge"] = [region_meta["on_edge"]]
                    region_feature.qualifiers["product"] = [region_meta["bgc_class"]]
                    region_feature.qualifiers["region_number"] = [region_meta["num"]]
                    region_feature.qualifiers["tool"] = [tool]

                    new_features.append(region_feature)

                    # only keep CDS and genes, which is relevant for BiG-SLiCE
                    for feature in gbk.features:
                        if feature.type in ["CDS", "gene"]:
                            # print(feature)
                            if not feature.qualifiers.get("translation", None):
                              # print((
                              #     "CDS {} is not translated "
                              #     "(don't have 'translation' qualifier)!"
                              #     " translating with BioPython "
                              #     "(default transl_table=11 if not filled).."
                              # ).format(feature.qualifiers.get("locus_tag", "n/a")))
                                feature.qualifiers["transl_table"] = [feature.qualifiers.get("transl_table", 11)]
                                feature.qualifiers["translation"] = [feature.translate(gbk, cds=False).seq]
                            new_features.append(feature)

                    gbk.features = new_features

        ############################################################################
        # Write GBKs 
        ############################################################################ 

        file_name = os.path.basename(input_gbk)
        output_gbk = os.path.join(output_dir, file_name)
        SeqIO.write(gbk, output_gbk, "gb")

# format_gbks(input_dir = "../tests/OceanDNA-a1001/bgc_annot/inter/deepbgc/gbks", \
#             input_tsv = "../tests/OceanDNA-a1001/tables/deepbgc_annot_metadata.tsv", \
#             output_dir = "../tests/output_test/bgc_annot/inter/deepbgc/test/tmp2", \
#             tool  = "deepbgc")

# format_gbks(input_dir = "../tests/output_test/bgc_annot/inter/gecco/test_redu2", \
#             input_tsv = "../tests/output_test/tables/gecco_annot_metadata.tsv", \
#             output_dir = "../tests/output_test/bgc_annot/inter/gecco/test_redu2/gbk", \
            # tool  = "gecco")

'''
Description
This function takes as an input GenkBank files containing one BGC sequence per file and and the `*_annot_metadata.tsv` table (generated with the corresponding `*_annot_pareser` function). The output consists of GenBank files following antiSMASH specifications.
The main modifications are including the “product”, “region_numnber”, and “contig_edge” features, which are obtained from the `*_annot_metadata.tsv` files. 
Another relevant modification is that the `LOCUS` field is modified to have as a name the contig_id (in the case of deepBGC and gecco).
The code was adapted from a script provided by BiG-SLICE, to generate antiSMASH-like BGC sequences.
'''