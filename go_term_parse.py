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
import seqscreen_parse_utils

def main():
    """
    Runs the go term parsing from the command line
    takes in a .tsv file and a list of go terms as input
    """

    # parse the inputted .tsv file
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input a .tsv file")
    parser.add_argument(
        "go_numbers", type=str, nargs="+",
        help="input go numbers separated by spaces")
    args = parser.parse_args()
    filename = pathlib.PurePath(args.input_file).stem
    output_dir = "outputs/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # remove rows where there are no bpocs, create a revised file
    dataframe = pd.read_csv(pathlib.Path(args.input_file), sep='\t', dtype=str)
    go_nums = args.go_numbers
    for go_num in go_nums:
        seqscreen_parse_utils.go_term_parse(dataframe, go_num, filename, output_dir)
        krona_input = os.path.join(output_dir, filename + f"_go_num_{go_num}_revised.tsv")
        seqscreen_parse_utils.krona_plot(krona_input)


main()
