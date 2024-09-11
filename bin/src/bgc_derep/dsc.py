###############################################################################
## Function dereplicate_shared_contigs
###############################################################################

def dereplicate_shared_contigs(metadata_df1: str = None, 
                               metadata_df2: str = None,
                               dereplicated_bgcs: str = None,
                               overlapped_bgcs: str = None,
                               partially_overlapped_bgcs: str = None,
                               non_overlapped_bgcs: str = None):
    
    for tool in dereplicated_bgcs:
        dereplicated_bgcs[tool] = dict()

    '''
    Get get contigs interesection  
    '''
    shared_contigs_list = list(set(metadata_df1['contig_id']).intersection(metadata_df2['contig_id']))

    '''
    Iterate through contits  
    '''
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
                    overlap = e1 - s1
                    overlapped_bgcs[tool1].update({id1: {'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                    if id1 in dereplicated_bgcs[tool1]:
                        del dereplicated_bgcs[tool1][id1] # id1 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 2nd case:
                # tool1          s1----e1
                # tool2      s2----------------e2

                if s1 > s2 and e1 < e2:
                    case = 2
                    overlap = e1 - s1 +1
                    overlapped_bgcs[tool1].update({id1: {'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2: {'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file': f2, 'row': tool2_row}})
                    if id1 in dereplicated_bgcs[tool1]:
                        del dereplicated_bgcs[tool1][id1] # id1 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 3rd case:
                # tool1                 s1----e1
                # tool2       s2--------------e2

                if s1 > s2 and e1 == e2:
                    case = 3
                    overlap = e1 - s1 +1
                    overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file':f2, 'row': tool2_row}})
                    if id1 in dereplicated_bgcs[tool1]:
                        del dereplicated_bgcs[tool1][id1] # id1 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 4th case:
                # tool1             s1--------------e1
                # tool2      s2---------------e2

                if s1 > s2 and e1 > e2 and e2 > s1:
                    case = 4
                    overlap = e2 - s1 +1
                    if l1 >= l2:
                        lref = l1
                        dereplicated_bgcs[tool1].update({id1 :{'file':f1,'row':tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2] 
                            # id2 could exists if it was the longer seq in cases 1st,2nd, 10th, 11th, 6th or 12th, in a previous iteration
                    else:
                        lref = l2
                        dereplicated_bgcs[tool2].update({id2: {'file':f2,'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1] 
                            # id1 could exists if it was the longer seq in cases 8th, 9th, 10th, 11th, 6th or 12th in a previous iteration

                    if overlap/lref <= overlap_thres:
                        partially_overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                        partially_overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    else:
                        overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                        overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})

                # 5th case:
                # tool1                       s1--------------e1
                # tool2      s2---------------e2

                if s1 > s2 and e1 > e2 and e2 == s1:
                    case = 5
                    overlap = 1
                    partially_overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2,'file':f1, 'row': tool1_row}})
                    partially_overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})

                    if l1 >= l2:
                        dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2] 
                            # id2 could exists if it was the longer seq in cases 1st, 2dn, 10th,11th, 6th or 12th in a previous iteration
                    else:
                        dereplicated_bgcs[tool2].update({id2: {'file':f2, 'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1] 
                            # id1 could exists if it was the longer seq in cases 8th, 9th,10th, 11th, 6th or 12th in a previous iteration

                # 6th case:
                # tool1                            s1--------------e1
                # tool2      s2---------------e2

                if s1 > s2 and e1 > e2 and e2 < s1:
                    case = 6
                    overlap = None
                    non_overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    non_overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file':f2, 'row': tool2_row}})

                # 7th case:
                # tool1      s1---------------e1
                # tool2      s2------e2

                if s1 == s2 and e1 > e2:
                    case = 7
                    overlap = e2 - s2 +1
                    overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1,'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2,'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                    if id2 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id2] # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 8th case:
                # tool1    s1---------------e1
                # tool2         s2----e2

                if s1 < s2 and e1 > e2:
                    case = 8
                    overlap = e2 - s2 +1
                    overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                    if id2 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id2] # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 9th case:
                # tool1      s1---------------e1
                # tool2                s2-----e2

                if s1 < s2 and e1 == e2:
                    case = 9
                    overlap = e2 - s2 +1
                    overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                    if id2 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id2] # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

                # 10th case:
                # tool1      s1--------------e1
                # tool2               s2-------------e2

                if s1 < s2 and e1 < e2 and s2 < e1:
                    case = 10
                    overlap = e1 - s2 +1

                    if l1 >= l2:
                        lref = l1
                        dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2] 
                            # id2 could exists if it was the longer seq in cases 2nd, 3rd, 4th, 5th, 6th or 12th in a previous iteration
                    else:
                        lref = l2
                        dereplicated_bgcs[tool2].update({id2: {'file':f2, 'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1] 
                            # id1 could exists if it was the longer seq in cases 7th, 8th, 4th, 5th, 6th or 12th in a previous iteration

                    if overlap/lref <= overlap_thres:
                        partially_overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                        partially_overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    else:
                        overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                        overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})

                # 11th case:
                # tool1     s1--------------e1
                # tool2                     s2---------------e2

                if s1 < s2 and e1 < e2 and e1 == s2:
                    case = 11
                    overlap = 1
                    partially_overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    partially_overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})

                    if l1 >= l2:
                        dereplicated_bgcs[tool1].update({id1 :{'file':f1, 'row': tool1_row}})
                        if id2 in dereplicated_bgcs[tool2]:
                            del dereplicated_bgcs[tool2][id2] 
                            # id2 could exists if it was the longer seq in cases 2nd,3rd, 4th, 5th, 6th or 12th in a previous iteration
                    else:
                        dereplicated_bgcs[tool2].update({id2: {'file':f2, 'row': tool2_row}})
                        if id1 in dereplicated_bgcs[tool1]:
                            del dereplicated_bgcs[tool1][id1] 
                            # id2 could exists if it was the longer seq in cases 7th, 8th, 4th, 5th, 6th or 12th in a previous iteration

                # 12th case:
                # tool1     s1--------------e1
                # tool2                           s2---------------e2

                if s1 < s2 and e1 < e2 and e1 < s2:
                    case = 12
                    overlap = None
                    non_overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    non_overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                    dereplicated_bgcs[tool2].update({id2: {'file':f2, 'row': tool2_row}})

                # 13th case:
                # tool1      s1--------------e1
                # tool2      s2--------------e2

                if s1 == s2 and e1 == e2:
                    case = 13
                    overlap = e1 - s1 +1
                    overlapped_bgcs[tool1].update({id1 :{'length':l1, 'overlap':overlap, 'case':case, 'match':id2, 'file':f1, 'row': tool1_row}})
                    overlapped_bgcs[tool2].update({id2 :{'length':l2, 'overlap':overlap, 'case':case, 'match':id1, 'file':f2, 'row': tool2_row}})
                    dereplicated_bgcs[tool1].update({id1: {'file':f1, 'row': tool1_row}})
                    if id1 in dereplicated_bgcs[tool2]:
                        del dereplicated_bgcs[tool2][id1] # id2 could exists if it was the longer seq in cases 6th or 12th in a previous iteration

    output_dict = {'dereplicated_bgcs': dereplicated_bgcs, 
                   'overlapped_bgcs': overlapped_bgcs, 
                   'partially_overlapped_bgcs': partially_overlapped_bgcs, 
                   'non_overlapped_bgcs': non_overlapped_bgcs}

    return(output_dict)

