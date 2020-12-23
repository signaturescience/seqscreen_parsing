#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
go_term_parse.py

Script to extract information regarding GO terms in sequences.
"""
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
    print (go_nums)
    for go_num in go_nums:
        if go_num[:3] == "GO:":
            go_num = go_num[3:]
        print(go_num,"next")
        # Must add logic to NOT continue if no go_terms are found.  (logic goes in go_term_parse, add catch here)
        seqscreen.go_term_parse(dataframe, go_num, filename, output_dir)
#        krona_input = os.path.join(output_dir, filename + f"_go_num_{go_num}_revised.tsv")
#        status = seqscreen.krona_plot(krona_input)
#        if status != 0:
#            print (status, go_num, "Problem")


main()
