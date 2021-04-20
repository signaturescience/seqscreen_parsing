#!/usr/bin/env python

####
# Created by Matthew Scholz for Signature Science
# 12/21/2020
####

import goatools
from goatools.obo_parser import GODag
from goatools.godag.go_tasks import get_go2parents
import argparse
import os
import pathlib
import pandas as pd
import argparse
import sys
import pathlib
import pandas as pd
import seqscreen_parse_utils as seqscreen
import re
import numpy as np
import subprocess
import timeit

def main():
    """
    Runs the go term parsing from the command line
    takes in a .tsv file and a list of go terms as input
    """
    # parse the inputted .tsv file
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input a .tsv file")
    #    parser.add_argument("go_numbers", type=list, nargs="?", help="input go numbers separated by spaces")
    parser.add_argument("-g", "--go_file", required=False, type=str, help="File containing GO terms to be parsed")
    parser.add_argument("-G", "--go_terms", required=False, type=str, action='append', help="Add individual GO terms can be repeated multiple times (e.g. -G GO:000112 -G GO:000013)" )
    parser.add_argument("--prefix", "-p", type=str, help="prefix for output files, default is seqscreen_GO", default="seqscreen_GO")
    parser.add_argument("--out","-o", type=str, required=False, help="output directory, default = output", default="output")
    parser.add_argument("--krona", "-k", action='store_true', help="generate krona plots per GO term")
    parser.add_argument("--force", "-f", help="Force run if files exist", action='store_true')
    parser.add_argument("--quiet", "-q", help="run quietly, minimal printing to STDOUT", action='store_true')
    args = parser.parse_args()
    if not args.force:
        check_file = "".join((args.out,"/",args.prefix,".full_report.csv"))
        print(check_file)
        if os.path.exists(check_file):
            print("file: ",args.out,args.prefix,".full_report.csv"+" exists, use --force/-f to force overwrite", sep="")
            exit(1)
    filename = pathlib.PurePath(args.input_file).stem
    #check for output directory, create if not exists
    output_dir = args.out
    seqscreen.create_output_directory(output_dir)
    #all GO terms need to be concatenated
    go_nums = []
    if args.go_terms is not None:
        for go_term in args.go_terms:
            go_nums.append(go_term)
    if args.go_file is not None and os.path.exists(args.go_file):
        with open(args.go_file, 'r') as go_reader:
            for line in go_reader:
                go_nums.append(line.strip())
    if args.quiet:
        f=open(os.devnull, "w")
        sys.stdout = f
    if len(go_nums) == 0:
        print(args.go_numbers,"and", args.go_file, "are blank, please either define go terms on the command line, or provide a file")
        parser.print_help()
        exit(1)
    # remove rows where there are no bpocs, create a revised file
    dataframe = pd.read_csv(pathlib.Path(args.input_file), sep='\t', dtype=str)
    #fix filename to variable
    godag = GODag("go-basic.obo", optional_attrs={'relationship'})
    #parsed_dataframe = seqscreen.parse_GO_terms(godag,dataframe,go_nums)
    # FOR TROUBLESHOOTING
#parsed_dataframe = parse_GO_terms(godag, dataframe, go_nums)
    parsed_dataframe = seqscreen.parse_GO_terms(godag, dataframe, go_nums)
    #parsed_dataframe = numpy_parse_GO_terms(godag, dataframe, go_nums)
    grouped_list = seqscreen.collapse_GO_results(parsed_dataframe,['GO_term','taxid'], 'taxid')
    taxID_list = seqscreen.collapse_GO_results(parsed_dataframe, ['taxid'], 'taxid')
    uniprot_list = seqscreen.collapse_GO_results(parsed_dataframe,['uniprot'], 'uniprot')
    uniprot_list.to_csv("{OUT}/{PRE}.uniprot.csv".format(PRE=args.prefix, OUT=args.out))
    grouped_list.to_csv("{OUT}/{PRE}.grouped.csv".format(PRE=args.prefix, OUT=args.out))
    taxID_list.to_csv("{OUT}/{PRE}.taxids.csv".format(PRE=args.prefix, OUT=args.out))
    parsed_dataframe.to_csv("{OUT}/{PRE}.full_report.csv".format(PRE=args.prefix, OUT=args.out))
    # iterate through GO_Terms to make krona plot for each:
    #krona_format=['query','taxid']
    if args.krona:
        for go in go_nums:
            sliced=parsed_dataframe.loc[parsed_dataframe['GO_term'] == go]
            name = godag[go].name
            name = re.sub('\W+',"_",name)
            if len(sliced) > 0:
                sliced = sliced.groupby(['query','taxid'])[['taxid']].count()
                go_filename = go.replace(":","-")
                file_root = "{OUT}/{PRE}.{GO}.{NAME}.krona".format(NAME=name, PRE=args.prefix, GO=go_filename, OUT=args.out)
                seqscreen.krona_from_slice(sliced,file_root, name)





def numpy_parse_GO_terms(godag, dataframe, go_nums):
    return_df = pd.DataFrame(columns=['GO_term', 'query', 'organism', 'associated_GO_terms', 'multi_taxids_confidence',
                                      'taxid', 'gene_name', 'uniprot', 'uniprot evalue'], dtype='str')
    """
    This structure is horribly inefficient.  Needs to be refactored to do the following:
    generate expanded dataframe of ALL queries and ALL GO terms
    Filter to only include those rows that contain GO term
    check query/GO combo for any GO term children, place in associated_GO_terms
    """
    rows_raw = len(dataframe.index)
    dataframe2 = dataframe
    (dataframe, num_rows) = seqscreen.remove_blanks_from_dataframe(dataframe, 'go')
    print(num_rows, "of", rows_raw, "with annotated GO terms kept")
    dataframe['go_id_confidence'] = dataframe['go_id_confidence'].str.split(";")
    expanded_dataframe = dataframe.explode('go_id_confidence')
    expanded_dataframe['go'] = expanded_dataframe['go_id_confidence'].str.replace("\[.*", "", regex=True)
    godag_keys = godag.keys()
    total_iter = 0
    go_rows = np.asanyarray(expanded_dataframe['go'])
    for go in go_nums:
        print(go)
        string_iter = 0
        #slice_go = expanded_dataframe.loc[expanded_dataframe['go'] == go]
        funny = list(np.flatnonzero(go_rows == go))
        print(len(go_rows))
        print(max(funny))
        print(len(expanded_dataframe['go']))
        slice_go = expanded_dataframe.iloc[funny]
        go_family = [go]
        if go in godag_keys:
            for value in list(godag[go].get_all_children()):
                go_family.append(value)
        if len(slice_go.index) > 0:
            print("0", "/", len(slice_go.index))
            queries = list(set(slice_go['query']))
            string_iter = 0
            orig_slice = np.asanyarray(dataframe2['query'])
            df2_index = np.in1d(orig_slice, queries)
            orig_slice = np.asanyarray(dataframe['query'])    
            df_index = np.intersect1d(orig_slice, queries)
            df2_subset = dataframe2.iloc[df2_index]
            orig_slice = np.asanyarray(expanded_dataframe['query'])
            edf_slice = np.intersect1d(orig_slice, queries)
            edf_subset = expanded_dataframe.iloc[edf_slice]
            for query in queries:
                total_iter += 1
                string_iter += 1
                if string_iter % 1000 == 0:
                    print(string_iter, "/", len(queries))
                fl = edf_subset['go']['query' == query and 'go' in go_family]
                query_list=np.intersect1d(fl,go_family)
                idx = slice2.index[0]
                return_df.loc[total_iter] = {'GO_term': go, 'query': query,
                                             'organism': slice2['organism'][idx],
                                             'associated_GO_terms': ";".join(query_list),
                                             'multi_taxids_confidence': slice2['multi_taxids_confidence'],
                                             'taxid': slice2['taxid'][idx],
                                             'gene_name': slice2['gene_name'][idx],
                                             'uniprot': slice2['uniprot'][idx],
                                             'uniprot evalue': slice2['uniprot evalue'][idx]}  # dicer
        print(string_iter, "/", string_iter, ":", go, "Complete")
    return_df.to_csv("test.csv")
    return (return_df)



if __name__ == "__main__":
    # execute only if run as a script
    sys.exit(main())

