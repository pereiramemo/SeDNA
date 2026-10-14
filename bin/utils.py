###############################################################################
# 1. Set env
###############################################################################

# Only the standard library is imported at module level. utils.py is shared by
# every module image, and those images do not share a dependency set: module 2
# (deepBGC) runs Python 3.7, and get_bgc_coverage.py runs without biopython.
# pandas and Bio are therefore imported inside the functions that need them,
# and the syntax here must stay Python 3.7 compatible (typing.List, no walrus,
# no list[str] generics).

import gzip
import os
import re
import shutil
import subprocess
import sys
from typing import List, Optional

###############################################################################
# 2. Generic utility functions
###############################################################################

###############################################################################
# 2.1 Run a command and check return code; optionally redirect stdout to file
###############################################################################

def run(cmd: List[str], stdout_path: Optional[str] = None) -> None:
    """Run a command, stream output or redirect to file, and fail on non-zero return code."""
    try:
        if stdout_path:
            with open(stdout_path, "w") as out_f:
                proc = subprocess.run(cmd, stdout=out_f, check=False)
        else:
            proc = subprocess.run(cmd, check=False)
    except FileNotFoundError as exc:
        print(f"Command not found: {cmd[0]} ({exc})", file=sys.stderr)
        sys.exit(1)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd)

###############################################################################
# 2.2 Check that required tools are available on PATH
###############################################################################

def check_tools(tools: List[str]) -> None:
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        print(f"Missing tools: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.3 Check a file exists; exit with error if not
###############################################################################

def check_file(path: str, label: str) -> None:
    if not os.path.isfile(path):
        print(f"{label} is not a real file", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.4 Check that a FASTA file has at least one record
###############################################################################

def check_fasta(path: str) -> None:
    from Bio import SeqIO
    try:
        opener = gzip.open if is_gzipped(path) else open
        with opener(path, "rt") as handle:
            fasta_check = any(SeqIO.parse(handle, "fasta"))
        if not fasta_check:
            raise ValueError(f"{path} contains no valid FASTA records.")
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.5 Create an output directory, honouring --overwrite
###############################################################################

def prepare_output_dir(path: str, overwrite: bool) -> None:
    """Create path; if it exists, remove it first when overwrite is set, else exit."""
    if os.path.exists(path):
        if not overwrite:
            print(f"Output dir {path} already exists. Use --overwrite to overwrite",
                  file=sys.stderr)
            sys.exit(1)
        try:
            shutil.rmtree(path)
        except OSError as e:
            print(f"Error removing {path}: {e}", file=sys.stderr)
            sys.exit(1)
    try:
        os.makedirs(path)
    except Exception as e:
        print(f"Error creating directory {path}: {e}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2.6 Detect gzip compression by magic bytes
###############################################################################

def is_gzipped(path: str) -> bool:
    """True if path starts with the gzip magic number, whatever it is named."""
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"

###############################################################################
# 3. BGC annotation utilities (formerly bin/src/bgc_annot/utilities.py)
###############################################################################

###############################################################################
# 3.1 find_files
###############################################################################

def find_files(input_dir: str = None, pattern: str = r".*\.region\d+.gbk") -> list:
    '''
    Given a directory and pattern, recursively find all the files that match
    the pattern (a regular expression), and return a list with their paths.
    '''
    matching_files = []
    for root, dirs, files in os.walk(input_dir, followlinks=True):
        for file in files:
            if re.search(pattern, file):
                matching_files.append(os.path.join(root, file))
    return matching_files

###############################################################################
# 3.2 find_dirs
###############################################################################

def find_dirs(input_dir: str = None, pattern: str = r".*\.region\d+.gbk") -> list:
    '''
    Given a directory and pattern, recursively find all the dirs that match
    the pattern (a regular expression), and return a list with their paths.
    '''
    matching_dirs = []
    for root, dirs, files in os.walk(input_dir, followlinks=True):
        for dir in dirs:
            if re.search(pattern, dir):
                matching_dirs.append(os.path.join(root, dir))
    return matching_dirs

###############################################################################
# 3.3 get_features
###############################################################################

def get_features(record=None, feature_name: str = None, qualifier_name: str = None) -> str:
    '''
    Extract the qualifier `qualifier_name` within `feature_name` from a GBK
    record parsed with Bio.SeqIO.
    '''
    for feature in record.features:
        if feature.type == feature_name:
            return feature.qualifiers[qualifier_name][0]

###############################################################################
# 3.4 get_feature_location
###############################################################################

def get_feature_location(record=None, feature_name: str = None) -> list:
    '''
    Return the [start, end] coordinates of `feature_name` in a GBK record.
    '''
    for feature in record.features:
        if feature.type == feature_name:
            return [int(feature.location.start), int(feature.location.end)]

###############################################################################
# 3.5 get_seqs_length
###############################################################################

def get_seqs_length(input_fasta: str = None) -> dict:
    '''
    Return a dictionary with the sequence id as key and the sequence length
    as value, for a (optionally gzipped) FASTA file.
    '''
    from Bio import SeqIO
    sequences_lengths = {}
    opener = gzip.open if is_gzipped(input_fasta) else open
    with opener(input_fasta, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences_lengths[record.id] = len(record.seq)
    return sequences_lengths

###############################################################################
# 3.6 gbk_splitter
###############################################################################

def gbk_splitter(input_file: str = None, output_dir: str = None,
                 sample_name: str = None):
    '''
    Split a single GBK file containing multiple entries (deepBGC output) into
    one GBK per entry, named following the antiSMASH convention
    (<contig_id>.region<NNN>.gbk), and saved in output_dir.
    '''
    from Bio import SeqIO

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
        pattern1 = re.compile('_[0-9]+-[0-9]+$')
        pattern2 = re.compile(r'\.[0-9]$')
        if re.search(pattern2, acc):  # This is done to remove potential numbers after the coordinates. DeepBGC is inconsistent with this format.
            acc = pattern2.sub('', acc)
        contig_id = pattern1.sub('', acc)
        # create dict
        if contig_id not in bgc_contig_id_counter:
            bgc_contig_id_counter[contig_id] = 1
        else:
            bgc_contig_id_counter[contig_id] += 1
        bgc_region_names[acc] = "{0}.region{1:03d}.gbk".format(contig_id,
                                                               bgc_contig_id_counter[contig_id])

        output_file = os.path.join(output_dir, bgc_region_names[acc])
        SeqIO.write(gbk, output_file, "gb")

###############################################################################
# 3.7 gbk_filter
###############################################################################

def gbk_filter(input_dir: str = None, sample_name: str = None,
               deepbgc_score_thres: float = 0.75, deepbgc_cds_count_thres: int = 2,
               output_dir: str = None):
    '''
    Filter out deepBGC BGC sequences having a low score or a low number of
    CDSs. Input files must follow the antiSMASH naming (.*region.*.gbk).
    Low confidence BGCs are moved from input_dir to output_dir.
    '''
    from Bio import SeqIO

    def count_cds(record):
        cds_count = 0
        for feature in record.features:
            if feature.type == "CDS":
                cds_count += 1
        return int(cds_count)

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    try:
        os.makedirs(output_dir)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

    paths2gbk = find_files(input_dir=input_dir, pattern=".*region.*.gbk")
    for input_gbk in paths2gbk:
        with open(input_gbk, "r") as gbk_handle:
            record = SeqIO.read(gbk_handle, "genbank")

        # get deepbgc score
        deepbgc_score = get_features(record=record,
                                     feature_name="cluster",
                                     qualifier_name="deepbgc_score")
        deepbgc_score = float(deepbgc_score)

        # count number of CDS
        cds_count = count_cds(record=record)

        # Remove low confidence BGCs
        if deepbgc_score < deepbgc_score_thres or cds_count < deepbgc_cds_count_thres:
            shutil.move(input_gbk, output_dir)

###############################################################################
# 3.8 antismash_annot_parser
###############################################################################

def antismash_annot_parser(input_dir: str = None, sample_name: str = None,
                           tool: str = "antismash", input_fasta: str = None) -> "pandas.DataFrame":
    '''
    Detect all antiSMASH GBK files (region[0-9]+.gbk) and parse each one into
    a DF with columns: acc, bgc_class, start, end, length, on_edge, contig_id,
    contig_length, file, sample_name and tool. The assembly FASTA is needed to
    get contig lengths. Since antiSMASH reuses the contig name as acc, rows are
    indexed by acc + a counter.
    '''
    import pandas as pd
    from Bio import SeqIO

    def gbk_parser(paths2gbk: list = None, fasta_lengths_dict: dict = None):
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
                bgc_class = get_features(record, "cand_cluster", "product")
                # start and end coordinates and length
                location = get_feature_location(record, "cand_cluster")  # these coords are zero based
                # The region GBK holds only the region's sequence, so feature
                # coords are relative to it; Orig. start (zero based) is the
                # region's start on the contig
                offset = int(record.annotations['structured_comment']['antiSMASH-Data']['Orig. start'])
                start = offset + int(location[0])
                end = offset + int(location[1])
                length = end - start  # 1 is not added given that coords are zero based
                # on edge
                on_edge = get_features(record, "cand_cluster", "contig_edge")
                # input gbk
                input_gbk = os.path.abspath(input_gbk)
                # add index
                acc_i = acc + "-" + str(i)
                i += 1
                # create output dict
                acc2x[acc_i] = {"acc": acc, "bgc_class": bgc_class,
                                "start": start, "end": end, "length": length,
                                "on_edge": on_edge, "contig_id": contig_id,
                                "contig_length": contig_length, 'file': input_gbk,
                                "sample_name": sample_name, 'tool': tool}
        return acc2x

    # find all antismash BGCs
    antismash_gbks = find_files(input_dir=input_dir, pattern=r'.region[0-9]+.gbk')
    if len(antismash_gbks) == 0:
        return

    # get sequences length
    fasta_lengths_dict = get_seqs_length(input_fasta)

    # get all antismash BGC ids and coords
    antismash_gbks_dict = gbk_parser(antismash_gbks, fasta_lengths_dict)

    # convert to df
    return pd.DataFrame.from_dict(antismash_gbks_dict, orient='index')

###############################################################################
# 3.9 deepbgc_annot_parser
###############################################################################

def deepbgc_annot_parser(input_dir: str = None, sample_name: str = None,
                         tool: str = "deepbgc", input_fasta: str = None) -> "pandas.DataFrame":
    '''
    Same output as antismash_annot_parser() for deepBGC GBKs. on_edge is
    derived from coordinates: start < 3 or end within 30 bp of the contig end.
    The BGC class is the product_class with the highest score.
    '''
    import pandas as pd
    from Bio import SeqIO

    def gbk_parser(paths2gbk: list = None, fasta_lengths_dict: dict = None):
        acc2x = dict()
        i = 1
        for input_gbk in paths2gbk:
            with open(input_gbk, "r") as gbk_handle:
                record = SeqIO.read(gbk_handle, "genbank")

                # acc
                acc = record.annotations['accessions'][0]

                # contig_id
                pattern1 = re.compile('_[0-9]+-[0-9]+$')
                pattern2 = re.compile(r'\.[0-9]$')
                if re.search(pattern2, acc):  # This is done to remove potential numbers after the coordinates. DeepBGC is inconsistent with this format.
                    acc = pattern2.sub('', acc)
                contig_id = pattern1.sub('', acc)

                # start and end coords, and length
                coords = get_feature_location(record, "cluster")  # these coords are zero based
                # The record holds only the BGC's sequence, so feature coords are
                # relative to it; the acc suffix _<start>-<end> (zero based) gives
                # the BGC's position on the contig
                offset = int(re.search('_([0-9]+)-[0-9]+$', acc).group(1))
                start = offset + int(coords[-2])
                end = offset + int(coords[-1])
                length = end - start  # 1 is not added given that coords are zero based

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
                acc2x[acc_i] = {"acc": acc, "bgc_class": bgc_class,
                                "start": start, "end": end, "length": length,
                                "on_edge": on_edge, "contig_id": contig_id,
                                'contig_length': contig_length, 'file': input_gbk,
                                "sample_name": sample_name, 'tool': tool}
        return acc2x

    # find all deepbgc BGCs
    deepbgc_gbks = find_files(input_dir=input_dir, pattern=r'.region[0-9]+.gbk')
    if len(deepbgc_gbks) == 0:
        return

    # get sequences length
    fasta_lengths_dict = get_seqs_length(input_fasta)

    # get all deepbgc BGC ids and coords
    deepbgc_gbks_dict = gbk_parser(deepbgc_gbks, fasta_lengths_dict)

    # convert to df
    return pd.DataFrame.from_dict(deepbgc_gbks_dict, orient='index')

###############################################################################
# 3.10 gecco_annot_parser
###############################################################################

def gecco_annot_parser(input_dir: str = None, sample_name: str = None,
                       tool: str = "gecco", input_fasta: str = None) -> "pandas.DataFrame":
    '''
    Same output as antismash_annot_parser() for GECCO, built from GECCO's own
    *.clusters.tsv rather than from GBKs. on_edge is derived as in
    deepbgc_annot_parser(). The `file` field points to the GBKs renamed to the
    antiSMASH convention in the sibling `gbks` folder, which must already exist.
    '''
    import pandas as pd

    def tbl_parser(path2tsv: str = None, fasta_lengths_dict: dict = None):
        acc2x = dict()
        i = 1
        df = pd.read_csv(path2tsv, sep='\t')
        for index, row in df.iterrows():
            # acc
            acc = row['cluster_id']
            # contig_id
            contig_id = row['sequence_id']
            # start and end coords, and length
            start = int(row['start']) - 1  # coords are converted to zero based
            end = int(row['end'])
            length = end - start  # 1 is not added given that coords are zero based
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
            acc_renamed = re.sub(r'_cluster_\d+', '.region{:03d}'.format(extracted_number), acc)
            input_gbk_dir = os.path.join(os.path.dirname(input_dir), "gbks")
            input_gbk_list = find_files(input_dir=input_gbk_dir, pattern=f'{acc_renamed}.gbk')
            input_gbk = input_gbk_list[0]
            input_gbk = os.path.abspath(input_gbk)
            # create output dict
            acc2x[acc_i] = {"acc": acc, "bgc_class": bgc_class,
                            "start": start, "end": end, "length": length,
                            "on_edge": on_edge, "contig_id": contig_id,
                            'contig_length': contig_length, 'file': input_gbk,
                            "sample_name": sample_name, 'tool': tool}
        return acc2x

    # find gecco cluster table
    gecco_tsv = find_files(input_dir=input_dir, pattern='.clusters.tsv')
    if len(gecco_tsv) == 0:
        return

    # get sequences length
    fasta_lengths_dict = get_seqs_length(input_fasta)

    # get all gecco BGC ids and coords
    gecco_tbl_dict = tbl_parser(gecco_tsv[0], fasta_lengths_dict)

    # convert to df
    return pd.DataFrame.from_dict(gecco_tbl_dict, orient='index')

###############################################################################
# 3.11 format_gbks
###############################################################################

def format_gbks(input_dir: str = None, input_tsv: str = None,
                tool: str = None, output_dir: str = None):
    '''
    Rewrite single-BGC GBK files as antiSMASH-like region GBKs, adding the
    "product", "region_number" and "contig_edge" qualifiers from the
    annot_metadata.tsv table, and setting LOCUS to the contig_id for deepBGC
    and GECCO. Adapted from a BiG-SLiCE script that generates antiSMASH-like
    BGC sequences.
    '''
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    # Create output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Parse coordinates file: start and end positions
    bgc_metadata = {}
    df = pd.read_csv(input_tsv, sep='\t')
    for idx, line in df.iterrows():

        acc = line['acc']
        start = int(line['start'])
        end = int(line['end'])
        bgc_class = line['bgc_class']
        on_edge = line['on_edge']

        if acc not in bgc_metadata:
            bgc_metadata[acc] = []

        bgc_metadata[acc].append({
            "start": start,
            "end": end,
            "bgc_class": bgc_class,
            "on_edge": on_edge,
            "num": idx
        })

    # Create GBK
    paths2gbk = find_files(input_dir=input_dir, pattern=r'.region[0-9]+.gbk')

    for input_gbk in paths2gbk:
        with open(input_gbk, "r") as gbk_handle:

            gbk = SeqIO.read(gbk_handle, "genbank")
            acc = gbk.annotations['accessions'][0]

            if tool == "deepbgc":

                pattern1 = re.compile('_[0-9]+-[0-9]+$')
                pattern2 = re.compile(r'\.[0-9]$')
                if re.search(pattern2, acc):  # This is done to remove potential numbers after the coordinates. DeepBGC is inconsistent with this format.
                    acc = pattern2.sub('', acc)
                contig_id = pattern1.sub('', acc)
                gbk.name = contig_id

            if tool == "gecco":
                contig_id = re.sub(r'_cluster_[0-9]+$', '', acc)
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
                            if not feature.qualifiers.get("translation", None):
                                feature.qualifiers["transl_table"] = [feature.qualifiers.get("transl_table", 11)]
                                feature.qualifiers["translation"] = [feature.translate(gbk, cds=False).seq]
                            new_features.append(feature)

                    gbk.features = new_features

        # Write GBKs
        file_name = os.path.basename(input_gbk)
        output_gbk = os.path.join(output_dir, file_name)
        SeqIO.write(gbk, output_gbk, "gb")

###############################################################################
# 3.12 extract_cds
###############################################################################

def extract_cds(input_gbk: str = None, cds_id: str = 'locus_tag', output_fasta: str = None):
    '''
    Extract the translated CDSs of a GenBank file into a multi-FASTA file,
    using the `cds_id` qualifier as the sequence identifier.
    '''
    from Bio import SeqIO

    cds_sequences = []
    with open(input_gbk, "r") as gb_handle:
        for record in SeqIO.parse(gb_handle, "genbank"):
            acc = record.annotations['accessions'][0]
            cds_n = 1
            for feature in record.features:
                if feature.type == "CDS":
                    # Extract CDS sequence and add it to the list
                    if feature.qualifiers.get(cds_id) is not None:
                        cds_id_instance = feature.qualifiers[cds_id][0]
                    else:
                        cds_id_instance = f'{cds_id}-notfound'

                    coords_start = feature.location.start
                    coords_end = feature.location.end
                    cds_header = f"{cds_id_instance}|{coords_start}-{coords_end}|{cds_n}"
                    cds_sequences.append((cds_header, feature.qualifiers['translation'][0]))
                    cds_n += 1

    # Write all CDS sequences to the output multi-FASTA file
    with open(output_fasta, "w") as fasta_handle:
        for gene_name, sequence in cds_sequences:
            fasta_handle.write(f">{acc}|{gene_name}\n")
            fasta_handle.write(f"{sequence}\n")

###############################################################################
# 3.13 resolve_domain_structure
###############################################################################

def resolve_domain_structure(path2names: list) -> None:
    '''
    Run cath-resolve-hits on <file>.domtblout for every base path in
    path2names, writing <file>_annotdoms_resolved.tsv.
    '''
    for file in path2names:
        domtblout = f'{file}.domtblout'
        output_resolved_tsv = f'{file}_annotdoms_resolved.tsv'

        command = [
            "cath-resolve-hits",
            "--input-format", "hmmer_domtblout",
            domtblout
        ]

        with open(output_resolved_tsv, "w") as outfile:
            result = subprocess.run(
                command,
                stdout=outfile,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )

        if result.returncode != 0:
            print("Error executing cath-resolve-hits", file=sys.stderr)
            print("Error message:\n", result.stderr, file=sys.stderr)
            sys.exit(1)

###############################################################################
# 4. BGC dereplication utilities (formerly bin/src/bgc_derep/)
###############################################################################

###############################################################################
# 4.1 make_bgc_id: <tool>__<sample_name>__<contig_id>__<start>__<end>
###############################################################################

def make_bgc_id(row_as_dict: dict, prefix_field: str = 'tool') -> str:
    bgc_id_str_list = [row_as_dict[prefix_field], row_as_dict['sample_name'],
                       row_as_dict['contig_id'], row_as_dict['start'],
                       row_as_dict['end']]
    return '__'.join([str(s) for s in bgc_id_str_list])

###############################################################################
# 4.2 dereplicate_shared_contigs
###############################################################################

def dereplicate_shared_contigs(metadata_df1=None,
                               metadata_df2=None,
                               dereplicated_bgcs: dict = None,
                               overlapped_bgcs: dict = None,
                               partially_overlapped_bgcs: dict = None,
                               non_overlapped_bgcs: dict = None,
                               overlap_thres: float = None):

    for tool in dereplicated_bgcs:
        dereplicated_bgcs[tool] = dict()

    # Get contigs intersection
    shared_contigs_list = list(set(metadata_df1['contig_id']).intersection(metadata_df2['contig_id']))

    # BGCs that lost a comparison (to a longer, or equal, overlapping BGC). The
    # disjoint cases 6 and 12 add both BGCs to dereplicated_bgcs, which used to
    # bring back a BGC already superseded by another pair of the same contig;
    # losers are therefore removed once all pairs have been compared.
    lost_bgcs = set()

    # Iterate through contigs
    for contig in shared_contigs_list:

        # subset data frames
        tool1_bgcs_i = metadata_df1["contig_id"] == contig
        tool1_bgcs_df = metadata_df1[tool1_bgcs_i]
        tool1 = metadata_df1['tool'].iloc[0]

        tool2_bgcs_i = metadata_df2["contig_id"] == contig
        tool2_bgcs_df = metadata_df2[tool2_bgcs_i]
        tool2 = metadata_df2['tool'].iloc[0]

        #########################################################
        # Determine BGC relative positions and select BGCs
        #########################################################

        for tool1_index, tool1_row in tool1_bgcs_df.iterrows():
            for tool2_index, tool2_row in tool2_bgcs_df.iterrows():

                s1 = tool1_row['start']
                e1 = tool1_row['end']
                s2 = tool2_row['start']
                e2 = tool2_row['end']
                l1 = tool1_row['length']
                l2 = tool2_row['length']
                id1 = tool1_index
                id2 = tool2_index
                f1 = tool1_row['file']
                f2 = tool2_row['file']

                # check if coordinates are ordered
                if s1 > e1:
                    print("tool1 BGC coordinates are not in order")
                    sys.exit()

                if s2 > e2:
                    print("tool2 BGC coordinates are not in order")
                    sys.exit()

                # 1st case:
                # tool1      s1----e1
                # tool2      s2---------------e2

                if s1 == s2 and e1 < e2:
                    case = 1
                    lost_bgcs.add((tool1, id1))
                    overlap = e1 - s1
                    overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                    if id1 in dereplicated_bgcs[tool1]:
                        del dereplicated_bgcs[tool1][id1]  # id1 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 2nd case:
                # tool1          s1----e1
                # tool2      s2----------------e2

                if s1 > s2 and e1 < e2:
                    case = 2
                    lost_bgcs.add((tool1, id1))
                    overlap = e1 - s1 + 1
                    overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                    if id1 in dereplicated_bgcs[tool1]:
                        del dereplicated_bgcs[tool1][id1]  # id1 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 3rd case:
                # tool1                 s1----e1
                # tool2       s2--------------e2

                if s1 > s2 and e1 == e2:
                    case = 3
                    lost_bgcs.add((tool1, id1))
                    overlap = e1 - s1 + 1
                    overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                    if id1 in dereplicated_bgcs[tool1]:
                        del dereplicated_bgcs[tool1][id1]  # id1 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 4th case:
                # tool1             s1--------------e1
                # tool2      s2---------------e2

                if s1 > s2 and e1 > e2 and e2 > s1:
                    case = 4
                    lost_bgcs.add((tool2, id2) if l1 >= l2 else (tool1, id1))
                    overlap = e2 - s1 + 1
                    if l1 >= l2:
                        lref = l1
                        dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2]
                            # id2 could exists if it was the longer seq in cases 1st,2nd, 10th, 11th, 6th or 12th, in a previous iteration
                    else:
                        lref = l2
                        dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1]
                            # id1 could exists if it was the longer seq in cases 8th, 9th, 10th, 11th, 6th or 12th in a previous iteration

                    if overlap / lref <= overlap_thres:
                        partially_overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                        partially_overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    else:
                        overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                        overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})

                # 5th case:
                # tool1                       s1--------------e1
                # tool2      s2---------------e2

                if s1 > s2 and e1 > e2 and e2 == s1:
                    case = 5
                    lost_bgcs.add((tool2, id2) if l1 >= l2 else (tool1, id1))
                    overlap = 1
                    partially_overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    partially_overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})

                    if l1 >= l2:
                        dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2]
                            # id2 could exists if it was the longer seq in cases 1st, 2dn, 10th,11th, 6th or 12th in a previous iteration
                    else:
                        dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1]
                            # id1 could exists if it was the longer seq in cases 8th, 9th,10th, 11th, 6th or 12th in a previous iteration

                # 6th case:
                # tool1                            s1--------------e1
                # tool2      s2---------------e2

                if s1 > s2 and e1 > e2 and e2 < s1:
                    case = 6
                    overlap = None
                    non_overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    non_overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})

                # 7th case:
                # tool1      s1---------------e1
                # tool2      s2------e2

                if s1 == s2 and e1 > e2:
                    case = 7
                    lost_bgcs.add((tool2, id2))
                    overlap = e2 - s2 + 1
                    overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                    if id2 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id2]  # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 8th case:
                # tool1    s1---------------e1
                # tool2         s2----e2

                if s1 < s2 and e1 > e2:
                    case = 8
                    lost_bgcs.add((tool2, id2))
                    overlap = e2 - s2 + 1
                    overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                    if id2 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id2]  # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 9th case:
                # tool1      s1---------------e1
                # tool2                s2-----e2

                if s1 < s2 and e1 == e2:
                    case = 9
                    lost_bgcs.add((tool2, id2))
                    overlap = e2 - s2 + 1
                    overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                    if id2 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id2]  # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 10th case:
                # tool1      s1--------------e1
                # tool2               s2-------------e2

                if s1 < s2 and e1 < e2 and s2 < e1:
                    case = 10
                    lost_bgcs.add((tool2, id2) if l1 >= l2 else (tool1, id1))
                    overlap = e1 - s2 + 1

                    if l1 >= l2:
                        lref = l1
                        dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2]
                            # id2 could exists if it was the longer seq in cases 2nd, 3rd, 4th, 5th, 6th or 12th in a previous iteration
                    else:
                        lref = l2
                        dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1]
                            # id1 could exists if it was the longer seq in cases 7th, 8th, 4th, 5th, 6th or 12th in a previous iteration

                    if overlap / lref <= overlap_thres:
                        partially_overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                        partially_overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    else:
                        overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                        overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})

                # 11th case:
                # tool1     s1--------------e1
                # tool2                     s2---------------e2

                if s1 < s2 and e1 < e2 and e1 == s2:
                    case = 11
                    lost_bgcs.add((tool2, id2) if l1 >= l2 else (tool1, id1))
                    overlap = 1
                    partially_overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    partially_overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})

                    if l1 >= l2:
                        dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2]
                            # id2 could exists if it was the longer seq in cases 2nd,3rd, 4th, 5th, 6th or 12th in a previous iteration
                    else:
                        dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1]
                            # id2 could exists if it was the longer seq in cases 7th, 8th, 4th, 5th, 6th or 12th in a previous iteration

                # 12th case:
                # tool1     s1--------------e1
                # tool2                           s2---------------e2

                if s1 < s2 and e1 < e2 and e1 < s2:
                    case = 12
                    overlap = None
                    non_overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    non_overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})

                # 13th case:
                # tool1      s1--------------e1
                # tool2      s2--------------e2

                if s1 == s2 and e1 == e2:
                    case = 13
                    lost_bgcs.add((tool2, id2))
                    overlap = e1 - s1 + 1
                    overlapped_bgcs[tool1].update({id1: {'length': l1, 'overlap': overlap, 'case': case, 'match': id2, 'file': f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length': l2, 'overlap': overlap, 'case': case, 'match': id1, 'file': f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file': f1, 'row': tool1_row}})
                    if id2 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id2]  # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

    for tool, bgc_id in lost_bgcs:
        dereplicated_bgcs[tool].pop(bgc_id, None)

    output_dict = {'dereplicated_bgcs': dereplicated_bgcs,
                   'overlapped_bgcs': overlapped_bgcs,
                   'partially_overlapped_bgcs': partially_overlapped_bgcs,
                   'non_overlapped_bgcs': non_overlapped_bgcs}

    return output_dict

###############################################################################
# 4.3 dereplicate_non_shared_contigs
###############################################################################

def dereplicate_non_shared_contigs(metadata_df1=None,
                                   metadata_df2=None,
                                   dereplicated_bgcs: dict = None,
                                   non_overlapped_bgcs: dict = None):

    # get contigs complement to intersection
    not_shared_contigs_list = list(set(metadata_df1['contig_id']).symmetric_difference(metadata_df2['contig_id']))
    not_shared_contigs_tool1_list = list(set(not_shared_contigs_list).intersection(metadata_df1['contig_id']))
    not_shared_contigs_tool2_list = list(set(not_shared_contigs_list).intersection(metadata_df2['contig_id']))

    tool1_bgcs_i = metadata_df1["contig_id"].isin(not_shared_contigs_tool1_list)
    tool1_bgcs_df = metadata_df1[tool1_bgcs_i]

    tool2_bgcs_i = metadata_df2["contig_id"].isin(not_shared_contigs_tool2_list)
    tool2_bgcs_df = metadata_df2[tool2_bgcs_i]

    # Populate the dereplicated_bgcs and non_overlapped_bgcs dictionaries
    for bgcs_df in [tool1_bgcs_df, tool2_bgcs_df]:
        for index, row in bgcs_df.iterrows():
            f = row['file']
            l = row['length']
            tool = row['tool']
            case = 14
            dereplicated_bgcs[tool].update({index: {'file': f, 'row': row}})
            non_overlapped_bgcs[tool].update({index: {'length': l, 'overlap': None, 'case': case,
                                                      'match': None, 'file': f, 'row': row}})

    output_dict = {'dereplicated_bgcs': dereplicated_bgcs,
                   'non_overlapped_bgcs': non_overlapped_bgcs}
    return output_dict

###############################################################################
# 4.4 sort_single_table
###############################################################################

def sort_single_table(metadata_df1=None,
                      dereplicated_bgcs: dict = None,
                      non_overlapped_bgcs: dict = None):

    # Populate the dereplicated_bgcs and non_overlapped_bgcs dictionaries
    for index, row in metadata_df1.iterrows():
        f = row['file']
        l = row['length']
        tool = row['tool']
        case = 14
        dereplicated_bgcs[tool].update({index: {'file': f, 'row': row}})
        non_overlapped_bgcs[tool].update({index: {'length': l, 'overlap': None, 'case': case,
                                                  'match': None, 'file': f, 'row': row}})

    output_dict = {'dereplicated_bgcs': dereplicated_bgcs,
                   'non_overlapped_bgcs': non_overlapped_bgcs}
    return output_dict

###############################################################################
# 4.5 copy_gbks
###############################################################################

def copy_gbks(input_dict: dict = None, output_dir: str = None) -> None:
    '''
    Copy each BGC's GBK into output_dir as <bgc_id>.gbk.

    These used to be symlinks into the per-tool output directories. Under
    Nextflow those live in upstream task work dirs, which are not mounted into
    this task's container and are gone after `nextflow clean`, so the catalog
    holds real files instead.
    '''
    for tool in input_dict:
        for bgc_id in input_dict[tool]:

            row_as_dict = input_dict[tool][bgc_id]['row'].to_dict()
            bgc_id_str = make_bgc_id(row_as_dict)
            source_path = input_dict[tool][bgc_id]['file']
            target_path = os.path.join(output_dir, f'{bgc_id_str}.gbk')

            if not os.path.exists(source_path):
                print(f'File path {source_path} in metadata table does not exist', file=sys.stderr)
                sys.exit(1)

            try:
                shutil.copyfile(source_path, target_path)
            except OSError as e:
                print(f"Error copying {source_path}: {e}", file=sys.stderr)
                sys.exit(1)

###############################################################################
# 4.6 create_df
###############################################################################

def create_df(input_dict: dict = None,
              output_dir: str = None,
              output_dir_gbks: str = None,
              df: bool = False,
              tsv: bool = True,
              id_prefix_field: str = 'tool'):
    '''
    Build the metadata table of a BGC dictionary, adding `bgc_id` and `link`.
    `link` is the GBK path relative to output_dir (gbks/<bgc_id>.gbk), so the
    published table stays valid wherever output_dir is moved. id_prefix_field
    names the column used as the bgc_id prefix (`tools` for core/extended).
    '''
    import pandas as pd

    output_df = pd.DataFrame()

    for tool in input_dict:
        for index_value in input_dict[tool]:

            row_as_dict = input_dict[tool][index_value]['row'].to_dict()
            input_as_row = input_dict[tool][index_value]['row'].to_frame().T

            # Create bgc_id
            bgc_id_str = make_bgc_id(row_as_dict, id_prefix_field)
            input_as_row['bgc_id'] = bgc_id_str

            # Create the link field, containing the path of the sequence
            link_path = os.path.join(output_dir_gbks, f'{bgc_id_str}.gbk')
            if not os.path.exists(link_path) and df is False:
                print(f"Error in {link_path}. The file does not exist", file=sys.stderr)
                sys.exit(1)
            input_as_row['link'] = os.path.relpath(link_path, output_dir)

            # Extend df
            output_df = pd.concat([output_df, input_as_row], ignore_index=True)

    # Optionally write dereplicated annot_metadata or output as a DF
    output_tsv = os.path.join(output_dir, 'annot_metadata.tsv')
    if df is True:
        return output_df
    if tsv is True:
        if output_df.shape[0] > 0:
            output_df.to_csv(output_tsv, sep='\t', index=False)

###############################################################################
# 4.7 map_synonyms
###############################################################################

def map_synonyms(synonyms_file: str = None, input_dict: dict = None) -> dict:
    '''
    Add `bgc_class_formatted` to every row, mapping each tool's BGC class onto
    the shared vocabulary in bgc_class_synonyms.tsv (column bgc_class_<tool>).
    '''
    import pandas as pd

    synonyms_df = pd.read_csv(synonyms_file, sep="\t")
    ref_column = 'bgc_class_formatted'
    synonyms_dicts = {}

    for tool in input_dict:
        for index_value in input_dict[tool]:

            row_as_df = input_dict[tool][index_value]['row']

            # The lookup column comes from the row's own `tool`, not the dict key:
            # from the second pairwise merge on, dereplicate_shared_contigs files
            # rows of the merged catalog under the key of its first row's tool,
            # so a deepBGC BGC can sit under 'antismash' (KeyError on 'Polyketide')
            row_tool = row_as_df.loc['tool']
            if row_tool not in synonyms_dicts:
                tool_column = f'bgc_class_{row_tool}'
                synonyms_dicts[row_tool] = synonyms_df.set_index(tool_column)[ref_column].to_dict()

            bgc_class_i = row_as_df.loc['bgc_class']
            row_as_df.loc['bgc_class_formatted'] = synonyms_dicts[row_tool][bgc_class_i]

            # update row value
            input_dict[tool][index_value]['row'] = row_as_df

    return input_dict

###############################################################################
# 4.8 find_bgc_groups: the predictions supporting each representative BGC
###############################################################################

def find_bgc_groups(preds: list, rep_ids: list) -> dict:
    '''
    Return {rep_id: [member bgc_id, ...]} with the representative first and
    the rest in input order.

    preds is every input prediction, of every tool, in input-table order, as
    dicts with bgc_id, contig_id, start and end. rep_ids are the bgc_ids kept
    by recursive_dereplication (the longest version). Two predictions overlap
    when they share or touch a coordinate (s1 <= e2 and s2 <= e1), which is
    what the 13 dereplication cases merge into one locus. A representative's
    members are the non-representative predictions reachable from it through
    overlaps without crossing another representative; reachability (not
    direct overlap) also catches a prediction that lost to an intermediate
    winner which in turn lost to the representative. Nothing here depends on
    which tools, or how many, produced the predictions.
    '''
    order = {p['bgc_id']: i for i, p in enumerate(preds)}
    rep_set = set(rep_ids)

    by_contig = {}
    for p in preds:
        by_contig.setdefault(p['contig_id'], []).append(p)

    neighbours = {p['bgc_id']: [] for p in preds}
    for contig_preds in by_contig.values():
        for i in range(len(contig_preds)):
            for j in range(i + 1, len(contig_preds)):
                a = contig_preds[i]
                b = contig_preds[j]
                if a['start'] <= b['end'] and b['start'] <= a['end']:
                    neighbours[a['bgc_id']].append(b['bgc_id'])
                    neighbours[b['bgc_id']].append(a['bgc_id'])

    groups = {}
    for rep_id in rep_ids:
        members = {rep_id}
        queue = [rep_id]
        while queue:
            current = queue.pop()
            for n in neighbours.get(current, []):
                if n in members or n in rep_set:
                    continue
                members.add(n)
                queue.append(n)
        others = sorted(members - {rep_id}, key=lambda m: order[m])
        groups[rep_id] = [rep_id] + others
    return groups

###############################################################################
# 4.9 most_supported_region: the core coordinates of a group
###############################################################################

def most_supported_region(intervals: list) -> tuple:
    '''
    Return (start, end, support) of the region covered by the most distinct
    tools, for intervals given as (start, end, tool) in zero-based half-open
    coordinates.

    When every tool overlaps every other, this is their intersection. When
    they only overlap in a chain, it is the region with the highest support;
    among equally supported regions the longest wins, then the leftmost. A
    tool with several intervals counts once. Valid for any number of tools.
    '''
    breakpoints = sorted(set([s for s, e, t in intervals] + [e for s, e, t in intervals]))

    # support of each elementary segment [b_i, b_i+1)
    segments = []
    for b1, b2 in zip(breakpoints[:-1], breakpoints[1:]):
        tools = set(t for s, e, t in intervals if s <= b1 and e >= b2)
        segments.append((b1, b2, len(tools)))

    if len(segments) == 0:
        # all intervals have zero length at the same point
        s, e, t = intervals[0]
        return s, e, len(set(t for s, e, t in intervals))

    max_support = max(support for b1, b2, support in segments)

    # merge adjacent max-support segments into runs
    runs = []
    for b1, b2, support in segments:
        if support != max_support:
            continue
        if runs and runs[-1][1] == b1:
            runs[-1][1] = b2
        else:
            runs.append([b1, b2])

    best = max(runs, key=lambda r: (r[1] - r[0], -r[0]))
    return best[0], best[1], max_support

###############################################################################
# 4.10 locate_record: a GBK record's origin on its contig
###############################################################################

def locate_record(contig_seq: str, record_seq: str, start: int) -> int:
    '''
    Return the zero-based contig position where record_seq begins.

    Found by sequence, not from tool-specific metadata, so that it works for
    any annotation tool: the record's own origin need not equal the BGC start
    (an antiSMASH region can be wider than its cand_cluster). Tries `start`
    first, then the occurrence nearest to it.
    '''
    contig_seq = contig_seq.upper()
    record_seq = record_seq.upper()
    if contig_seq[start:start + len(record_seq)] == record_seq:
        return start

    hits = []
    i = contig_seq.find(record_seq)
    while i != -1:
        hits.append(i)
        i = contig_seq.find(record_seq, i + 1)
    if len(hits) == 0:
        return None
    return min(hits, key=lambda h: abs(h - start))

###############################################################################
# 4.11 collect_cds: the CDSs of a group's members, on contig coordinates
###############################################################################

def collect_cds(member_gbks: list, contig_seq: str) -> list:
    '''
    Return the CDS features of the member GBKs, shifted to contig coordinates
    and deduplicated.

    member_gbks is a list of (gbk_path, start) in priority order (the
    representative first). Two CDSs on the same strand with the same stop
    coordinate are the same gene called by different gene finders; the copy
    of the higher-priority member is kept. A missing translation is computed,
    as format_gbks does.
    '''
    from Bio import SeqIO

    cds_features = []
    seen = set()
    for gbk_path, start in member_gbks:
        record = SeqIO.read(gbk_path, "genbank")
        origin = locate_record(contig_seq, str(record.seq), int(start))
        if origin is None:
            print(f"Error: the sequence of {gbk_path} was not found in its contig", file=sys.stderr)
            sys.exit(1)

        for feature in record.features:
            if feature.type != "CDS":
                continue
            if not feature.qualifiers.get("translation", None):
                feature.qualifiers["transl_table"] = [feature.qualifiers.get("transl_table", 11)]
                feature.qualifiers["translation"] = [str(feature.translate(record, cds=False).seq)]
            shifted = feature._shift(origin)
            strand = shifted.location.strand
            stop = int(shifted.location.start) if strand == -1 else int(shifted.location.end)
            if (strand, stop) in seen:
                continue
            seen.add((strand, stop))
            cds_features.append(shifted)

    return sorted(cds_features, key=lambda f: int(f.location.start))

###############################################################################
# 4.12 write_region_gbk: an antiSMASH-like region GBK for a contig interval
###############################################################################

def write_region_gbk(contig_id: str, contig_seq: str, start: int, end: int,
                     cds_features: list, bgc_class: str, on_edge: str,
                     tool: str, output_gbk: str) -> int:
    '''
    Write the contig interval [start, end) as an antiSMASH-like region GBK
    (same layout as format_gbks), holding the CDSs entirely inside it. Return
    the number of CDSs written.
    '''
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq(contig_seq[start:end]), id=contig_id, name=contig_id,
                       description=contig_id)
    record.annotations["molecule_type"] = "DNA"
    record.annotations["accessions"] = [contig_id]
    record.annotations["structured_comment"] = {
        "antiSMASH-Data": {
            "Version": "Simulated",
            "Orig. start": start,
            "Orig. end": end,
            "Note": "SeDNA dereplicated region built from several BGC predictions"}}

    features = [SeqFeature(FeatureLocation(0, end - start), type="protocluster")]
    region_feature = SeqFeature(FeatureLocation(0, end - start), type="region")
    region_feature.qualifiers["contig_edge"] = [on_edge]
    region_feature.qualifiers["product"] = [bgc_class]
    region_feature.qualifiers["tool"] = [tool]
    features.append(region_feature)

    n_cds = 0
    for feature in cds_features:
        if int(feature.location.start) >= start and int(feature.location.end) <= end:
            features.append(feature._shift(-start))
            n_cds += 1
    record.features = features

    with open(output_gbk, "w") as handle:
        from Bio import SeqIO
        SeqIO.write(record, handle, "genbank")
    return n_cds
