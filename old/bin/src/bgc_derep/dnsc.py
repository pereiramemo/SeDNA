###############################################################################
## Function dereplicate_non_shared_contigs
###############################################################################

def dereplicate_non_shared_contigs(metadata_df1: str = None, 
                                   metadata_df2: str = None,
                                   dereplicated_bgcs: str = None,
                                   non_overlapped_bgcs: str = None):


    # get contigs complement to interesection
    not_shared_contigs_list = list(set(metadata_df1['contig_id']).symmetric_difference(metadata_df2['contig_id']))
    not_shared_contigs_tool1_list = list(set(not_shared_contigs_list).intersection(metadata_df1['contig_id']))
    not_shared_contigs_tool2_list = list(set(not_shared_contigs_list).intersection(metadata_df2['contig_id']))

    tool1_bgcs_i = metadata_df1["contig_id"].isin(not_shared_contigs_tool1_list)
    tool1_bgcs_df = metadata_df1[tool1_bgcs_i]

    tool2_bgcs_i = metadata_df2["contig_id"].isin(not_shared_contigs_tool2_list)
    tool2_bgcs_df = metadata_df2[tool2_bgcs_i]
     
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

    for tool2_index, tool2_row in tool2_bgcs_df.iterrows():
        id2 = tool2_index
        f2 = tool2_row['file']
        l2 = tool2_row['length']
        tool2 = tool2_row['tool']
        case = 14
        dereplicated_bgcs[tool2].update({id2: {'file':f2, 'row': tool2_row}})
        non_overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':None, 'case':case, 
                                                 'match':None, 'file':f2, 'row': tool2_row}})

    output_dict = {'dereplicated_bgcs':dereplicated_bgcs, 
                   'non_overlapped_bgcs':non_overlapped_bgcs}                  
    return(output_dict)
