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
import pysam
import re

################################################################################
# 2. Parse parameters
################################################################################

parser = argparse.ArgumentParser(
    prog="import_annot_cds.py",
    description="Import annotated CDSs within BGC sequences."
)
parser.add_argument(
    "--annot_metadata", required=True,
    help="TSV file with BGC coordinates."
)
parser.add_argument(
    "--gff", required=True,
    help="GFF file containing gene features."
)
parser.add_argument(
    "--domtblout", required=True,
    help="domtblout file with domain annotations from hmmsearch."
)

parser.add_argument(
    "--overwrite", default = False, 
    action=argparse.BooleanOptionalAction, help="Overwrite output directory."
)

parser.add_argument(
    "--output_dir",
    help="Output directory to store filtered and resolved domain annotation files."
)
args = parser.parse_args()

annot_metadata = args.annot_metadata
gff = args.gff
domtblout = args.domtblout
overwrite = args.overwrite
output_dir = args.output_dir

################################################################################
# 3. Check if annot_metadata exists
################################################################################

if not os.path.exists(annot_metadata):
    print("annot_metadata file does not exists")
    sys.exit(0)

################################################################################
# 4. Create output dir
################################################################################

input_dir = os.path.dirname(annot_metadata)

if output_dir is None:
    output_dir = os.path.join(input_dir, "import_annot_cds")

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
# 5. Define functions
################################################################################

def check_and_index_gff(gff):
    """
    Check if the GFF file is indexed by Tabix (i.e. a .tbi file exists).
    If not, index the GFF file using Tabix.
    """

    gff_name, ext = os.path.splitext(gff)
    global sorted_gff
    global index_file
    sorted_gff = f'{gff_name}.sorted.gff.gz'
    index_file = f'{sorted_gff}.tbi'
    if not os.path.exists(index_file):
        print(f"Index file {index_file} not found. Sorting, compressing, and indexing GFF file...")

        #  sort
        sort_gff_command = ["bash",
                            "-c",
                            f"(grep '^#' {gff}; grep -v '^#' {gff} | sort -t \"$(printf '\\t')\" -k1,1 -k4,4n) | bgzip > {sorted_gff}"]

        sort_gff_result = subprocess.run(sort_gff_command,
                                         stdout=subprocess.PIPE, 
                                         stderr=subprocess.PIPE, 
                                         text=True)

        if sort_gff_result.returncode != 0:
            print(f"Error sorting GFF file:\n{sort_gff_result.stderr}")
            sys.exit(1)
        else:
            print("GFF file sorted and compressed successfully.")

        # index
        index_gff_command = ["tabix", "-p", "gff", sorted_gff]

        index_gff_result = subprocess.run(index_gff_command,
                                          stdout=subprocess.PIPE, 
                                          stderr=subprocess.PIPE, 
                                          text=True)

        if index_gff_result.returncode != 0:
            print(f"Error indexing GFF file:\n{index_gff_result.stderr}")
            sys.exit(1)
        else:
            print("GFF file indexed successfully.")
    else:
        print("GFF file already indexed.")


def extract_cds(gff, annot_metadata):
    """
    Extract CDS for specified regions from a Tabix-indexed file.

    Parameters:
        db (str): Path to the Tabix-indexed file.
        input_df (pandas.DataFrame): DataFrame where each row contains:
            - bgc_id
            - contig_id
            - start coordinate
            - end coordinate

    Returns:
        list: A list of lists, where each inner list contains:
              [bgc_id, contig_id, start, end, cds].
    """

    # open annot_metadata
    try:
        annot_metadata_df = pd.read_csv(annot_metadata, sep="\t")
    except Exception as e:
        print(f"Error reading metadata file {annot_metadata}: {e}")
        sys.exit(1)

    # Test if it has data
    if annot_metadata_df.shape[0] == 0:
        print(f"No BGCs found in {annot_metadata}")
        sys.exit(0)

    # check if the required columns are present
    required_columns = ["bgc_id", "contig_id", "start", "end"]
    for col in required_columns:
        if col not in annot_metadata_df.columns:
            print(f"Column '{col}' not found in metadata file {annot_metadata}.")
            sys.exit(1)

    # Create empty dict
    cds_ids_dict = {}

    # Open the Tabix file
    tabix_file = pysam.TabixFile(gff)

    # Loop over each row in the input DataFrame
    for row in annot_metadata_df.itertuples(index=False):

        bgc_id = row.bgc_id
        contig_id = row.contig_id
        start = row.start + 1 # convert coordinates as 1-based (as used in prodigal)
        end = row.end
        region = f'{contig_id}:{start}-{end}'

        # Fetch records from the Tabix file for the given region
        for record in tabix_file.fetch(region=region):

            record_list = record.split('\t')
            feature = record_list[2]
            if feature == "CDS":
                gene_id_list = record_list[8].split(';')
                cds_id = gene_id_list[0].replace('ID=','')

                # Initialized dictionary with empty list if necessary
                if bgc_id not in cds_ids_dict:
                    cds_ids_dict[bgc_id] = []

                # Append cds id
                cds_ids_dict[bgc_id].append(cds_id)

    tabix_file.close()
    return cds_ids_dict

def filter_domtblout_by_cds(domtblout, cds_ids_dict):
    """
    Filter the domtblout file (output of hmmsearch) to retain only rows for the given CDS gene IDs.
    Lines starting with '#' are passed through. For non-header lines, the first column (query id)
    is checked against the provided CDS gene IDs.

    Writes the filtered lines to pyhmmsearch_output.pfam_bgc_filt.domtblout.
    """

    # Initialize path2names
    path2names = []

    # Create a list of cds_ids to extract
    cds_ids_to_extract = {value for sublist in cds_ids_dict.values() for value in sublist}

    for bgc_id,cds_ids in cds_ids_dict.items():

        domtblout_bgc_filt_name = os.path.join(output_dir, bgc_id)
        domtblout_bgc_filt = os.path.join(output_dir, f'{domtblout_bgc_filt_name}.domtblout')
        path2names.append(domtblout_bgc_filt_name)

        try:
            with open(domtblout, "r") as infile, open(domtblout_bgc_filt, "w") as outfile:
                for line in infile:
                    if not line.startswith("#"):
                        parts = line.strip().split()
                        if parts[0] in cds_ids:
                            outfile.write(line)

        except Exception as e:
            print(f"Error filtering domtblout file: {e}")
            sys.exit(1)

    return(path2names)

def resolve_domain_structure(path2names):
    """
    Resolve domain structure for each file using cath-resolve-hits.

    This function constructs the necessary file names for each base file path in 'path2names'
    and runs the 'cath-resolve-hits' command. The command is defined as a list, and output
    redirection is handled by passing a file handle to the subprocess.

    Parameters:
        path2names (list): List of base file paths (without extensions) for each sample.
    
    Returns:
        None

    If the cath-resolve-hits command fails for any file (i.e. returns a non-zero exit code),
    an error message is printed and the program exits.
    """
    for file in path2names:
        domtblout = f'{file}.domtblout'
        output_resolved_tsv = f'{file}_annotdoms_resolved.tsv'

        # Build the command as a list
        command = [
            "cath-resolve-hits",
            "--input-format", "hmmer_domtblout",
            domtblout
        ]

        # Open the output file and run the command, redirecting stdout to the file
        with open(output_resolved_tsv, "w") as outfile:
            result = subprocess.run(
                command,
                stdout=outfile,
                stderr=subprocess.PIPE,
                text=True
            )

        if result.returncode != 0:
            print("Error executing cath-resolve-hits")
            print("Error message:\n", result.stderr)
            sys.exit(1)

################################################################################
# 6. Check input and output directories
################################################################################

if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
    except Exception as e:
        print(f"Error creating output directory {output_dir}: {e}")
        sys.exit(1)

################################################################################
# 7. Check and index the GFF file (if necessary)
################################################################################

check_and_index_gff(gff)

################################################################################
# 8. Extract CDS gene IDs from the GFF file using Tabix
################################################################################

if 'sorted_gff' in globals():
    gff = sorted_gff

cds_ids_dict = extract_cds(gff = gff, annot_metadata = annot_metadata) 

################################################################################
# 9. Filter domtblout file for the extracted CDS gene IDs
################################################################################

path2names = filter_domtblout_by_cds(domtblout = domtblout,
                                     cds_ids_dict = cds_ids_dict)

################################################################################
# 10. Run cath-resolve-hits
################################################################################

resolve_domain_structure(path2names = path2names)

################################################################################
# 11. End of execution
################################################################################

print("import_annot_cds module executed successfully")