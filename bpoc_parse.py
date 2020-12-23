#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bpoc_parse.py

Script to create report regarding BPoC data in sequences.
"""
import argparse
import os
import pathlib
import pandas as pd
import seqscreen_parse_utils
import sys


def main():
    """
    Runs the bpoc parsing from the command line
    takes in a .tsv file as input
    """
    #Command line Arguments (input file):
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input a .tsv file")
    args = parser.parse_args()
    filename = pathlib.PurePath(args.input_file).stem
    output_dir = "outputs/"

   #check for output_dir, create if not exists
    if not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir)
        except:
            print("could not create output directory:",output_dir)
            exit(1)

   # Check for file existing
    if not os.path.exists(args.input_file):
        print("Input file not found: "+args.input_file)
        parser.print_help()
        exit(1)
    # remove rows where there are no bpocs, create a revised file
    # parse the inputted .tsv file
    try:
        dataframe = pd.read_csv(pathlib.Path(args.input_file), sep='\t', dtype=str)
    except:
        print("error parsing file: "+args.input_file)
        print("input file should be .tsv file")
        parser.print_help()
        exit(1)

    #call bpoc_parse (I don't like the way this works)
    seqscreen_parse_utils.bpoc_parse(dataframe, filename, output_dir)
    krona_input = os.path.join(output_dir, filename + "_revised.tsv")
    #run krona_plot using revised.tsv file that is created in bpoc_parse (again, no bueno)
    seqscreen_parse_utils.krona_plot(krona_input)

if __name__ == '__main__':
    sys.exit(main())
