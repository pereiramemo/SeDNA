###############################################################################
## Function dereplicate_non_shared_contigs
###############################################################################

def sort_single_table(metadata_df1: str = None, 
                      dereplicated_bgcs: str = None,
                      non_overlapped_bgcs: str = None):


    tool1_bgcs_df = metadata_df1
     
    # Polulate the dereplicated_bgcs and non_overlapped_bgcs dictionaries
    for tool1_index, tool1_row in tool1_bgcs_df.iterrows():
        id1 = tool1_index
        f1 = tool1_row['file']
        l1 = tool1_row['length']
        tool1 = tool1_row['tool']
        case = 14
        dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
        non_overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':None, 'case':case, 
                                                 'match':None, 'file':f1, 'row': tool1_row}})

    output_dict = {'dereplicated_bgcs':dereplicated_bgcs, 
                   'non_overlapped_bgcs':non_overlapped_bgcs}                  
    return(output_dict)
