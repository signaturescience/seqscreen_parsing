# -*- coding: utf-8 -*-

import concurrent.futures
import functools
import pandas as pd
import numpy as np
import utils
import collections

def process_go_num(main_dict, expanded_dataframe, godag, hashed_query_dict, slice_go, go): 
    
    return_df_dict = collections.defaultdict(list)
    go_family = [go]
    if go in godag:
        # list concatenation
        go_family.extend(list(godag[go].get_all_children()))

    iteration = 0
    queries = []
    if len(slice_go) > 0:
        queries = list(set(slice_go['query']))
        for query in queries:
            if iteration % 1000 == 0:
                print(f'{go}: {iteration} / {len(queries)}')
            iteration += 1
            

            query_rows = hashed_query_dict[query]
            # pull slice of main_dictionary
            current_slice = utils.slice_dict(main_dict, hashed_query_dict[query])
            # slice2 = dataframe2.iloc[hashed_query_dict[query]]
            
            # fl = utils.slice_dict(expanded_dict, hashed_query_dict[query])
            fl = expanded_dataframe['go']['query' == query and 'go' in go_family]

            # find the intersection between the two lists           
            query_list = utils.intersection(fl, go_family)
            # query_list=np.intersect1d(fl,go_family)
            idx = 0
            # define new dict
            temp_dict = {'GO_term': go, 'query': query,
                        'organism': current_slice['organism'][idx],
                        'associated_GO_terms': ";".join(query_list),
                        'multi_taxids_confidence': current_slice['multi_taxids_confidence'],
                        'taxid': current_slice['taxid'][idx],
                        'gene_name': current_slice['gene_name'][idx],
                        'uniprot': current_slice['uniprot'][idx],
                        'uniprot evalue': current_slice['uniprot evalue'][idx]}
            
            # append key values to list of values within each dict key
            for k, v in temp_dict.items():
                return_df_dict[k].append(v)
    print(f'{go}: {iteration} / {len(queries)} COMPLETE')

    return return_df_dict
