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
import seqscreen_parse_utils as seqscreen
import argparse
import os
import pathlib
import pandas as pd
import seqscreen_parse_utils as seqscreen

def main():
    """
    Runs the go term parsing from the command line
    takes in a .tsv file and a list of go terms as input
    """

    # parse the inputted .tsv file
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input a .tsv file")
    parser.add_argument("go_numbers", type=list, nargs="?", help="input go numbers separated by spaces")
    parser.add_argument("-f", "--go_file", required=False, type=str, help="File containing GO terms to be parsed")
    parser.add_argument("--prefix", "-p", type=str, help="prefix for output files")
    args = parser.parse_args()
    filename = pathlib.PurePath(args.input_file).stem

    #check for output directory, create if not exists
    output_dir = "outputs/"
    seqscreen.output_directory(output_dir)

    #all GO terms need to be concatenated
    go_nums = []
    if args.go_numbers is not None:
        for go_term in args.go_numbers:
            go_nums.append(go_term)
    if args.go_file is not None:
        if os.path.exists(args.go_file):
            with open(args.go_file, 'r') as go_reader:
                for line in go_reader:
                    go_nums.append(line.strip())

    if len(go_nums) == 0:
        print(args.go_numbers,"and", args.go_file, "are blank, please either define go terms on the command line, or provide a file")
        parser.print_help()
        exit(1)

    # remove rows where there are no bpocs, create a revised file
    dataframe = pd.read_csv(pathlib.Path(args.input_file), sep='\t', dtype=str)

    godag = GODag("go-basic.obo", optional_attrs={'relationship'})
    parsed_dataframe = seqscreen.parse_GO_terms(godag,dataframe,go_nums)
    grouped_list = seqscreen.collapse_GO_results(parsed_dataframe)

    grouped_list.to_csv("{PRE}.grouped.csv".format(PRE=args.prefix))
    parsed_dataframe.to_csv("{PRE}.full_report.csv".format(PRE=args.prefix))





# optional_relationships = set()
# go2parents_isa = get_go2parents(godag, optional_relationships)
# print('{GO} parent: {P}'.format(
#     GO=GO_ID,
#     P=go2parents_isa[GO_ID]))
#
#
# optional_relationships = {'regulates', 'negatively_regulates', 'positively_regulates'}
# go2parents_reg = get_go2parents(godag, optional_relationships)
# print('{GO} parents: {P}'.format(
#     GO=GO_ID,
#     P=go2parents_reg[GO_ID]))
#
#
#
# from goatools.gosubdag.gosubdag import GoSubDag
#
# gosubdag_r0 = GoSubDag([GO_ID], godag, prt=None)
# print('{GO} ancestors: {P}'.format(
#     GO=GO_ID,
#     P=gosubdag_r0.rcntobj.go2ancestors[GO_ID]))

main()