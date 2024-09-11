###############################################################################
## Function create_links
###############################################################################
import os
import sys
import pandas as pd

def create_links(input_dict: str = None,
                 output_dir: str = None):

    for tool in input_dict:
        for bgc_id in input_dict[tool]:

            file_name = os.path.basename(input_dict[tool][bgc_id]['file'])
            row_as_dict = input_dict[tool][bgc_id]['row'].to_dict()
            tool_name = row_as_dict['tool']
            sample_name = row_as_dict['sample_name']
            contig_id = row_as_dict['contig_id']
            start_coord = row_as_dict['start']
            end_coord = row_as_dict['end']
            bgc_id_str_list = [tool_name, sample_name, contig_id, 
                              start_coord, end_coord]
            bgc_id_str_list = [str(s) for s in bgc_id_str_list]
            bgc_id_str = '__'.join(bgc_id_str_list)
            source_dir = input_dict[tool][bgc_id]['file']
            source_path = os.path.relpath(source_dir, output_dir)
            target_path = f'{output_dir}/{bgc_id_str}.gbk'

            if not os.path.exists(source_dir):
                print("File paths in metadata table do not exist")
                sys.exit()

            try:
                if os.path.exists(target_path):
                    os.unlink(target_path)
                os.symlink(source_path, target_path)
            except OSError as e:
                print(f"Error creating symbolic link: {e}")

       
    
###############################################################################
## Function create_dfs
###############################################################################
    
# extract dereplicated annot_metadata DF

def create_df(input_dict: str = None,
              output_dir: str = None,
              df = False,
              tsv = True):

    output_df = pd.DataFrame()

    for tool in input_dict:
        for index_value in input_dict[tool]:

            row_as_dict = input_dict[tool][index_value]['row'].to_dict()
            input_as_row = input_dict[tool][index_value]['row'].to_frame().T        

            # Create bgc_id
            tool_name = row_as_dict['tool']
            sample_name = row_as_dict['sample_name']
            contig_id = row_as_dict['contig_id']
            start_coord = row_as_dict['start']
            end_coord = row_as_dict['end']
            bgc_id_str_list = [tool_name, sample_name, contig_id, 
                               start_coord, end_coord]
            bgc_id_str_list = [str(s) for s in bgc_id_str_list]
            bgc_id_str = '__'.join(bgc_id_str_list)
            input_as_row['bgc_id'] = bgc_id_str

            # Crate the link filed, contaitnig the link of the sequence
            output_dir_path = os.path.abspath(output_dir)                
            link_name = f'{output_dir_path}/{bgc_id_str}.gbk'
            if not os.path.exists(link_name) and df == False:
                print(f"Error in {link_name}. The file does not exists")
                sys.exit()
            input_as_row['link'] = link_name

            # Extend df 
            output_df = pd.concat([output_df, input_as_row], ignore_index=True)

    # write dereplicated annot_metadata DF
    output_tsv = f'{output_dir}/annot_metadata.tsv'
    if df == True:
        return(output_df)
    if tsv == True: 
        if output_df.shape[0] > 0:
            output_df.to_csv(output_tsv, sep='\t', index=False)
