#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
taxid_parse.py

Script to extract information regarding taxon abundances in sequences.
"""
import argparse
import os
import pandas as pd
from seqscreen_parse_utils import *

def main():
    """
    Runs taxid parsing from command line.
    Takes in .tsv file as input.

    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("function",
                        nargs="?",
                        choices=['parse_conf', 'thresh_tied', 'all_tied',
                                 'assume_human', "count_taxid"],)
    args, sub_args = parser.parse_known_args()

    output_dir = "outputs/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if args.function == "parse_conf":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        parser.add_argument("confidence", type=float)
        args = parser.parse_args(sub_args)
        dframe = pd.read_csv(args.input, sep="\t")
        filename = (args.input.split("/")[-1]).split(".")[0]
        destname = f"{output_dir}{filename}_output.txt"
        parse_funcs(dframe, destname, sort_conf, args.confidence)
        make_krona(destname)

    elif args.function == "thresh_tied":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        parser.add_argument("threshold", type=int)
        args = parser.parse_args(sub_args)
        dframe = pd.read_csv(args.input, sep="\t")
        filename = (args.input.split("/")[-1]).split(".")[0]
        destname = f"{output_dir}{filename}_output.txt"
        parse_funcs(dframe, destname, sort_tied, args.threshold)
        make_krona(destname)

    elif args.function == "all_tied":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        args = parser.parse_args(sub_args)
        dframe = pd.read_csv(args.input, sep="\t")
        filename = (args.input.split("/")[-1]).split(".")[0]
        destname = f"{output_dir}{filename}_output.txt"
        parse_funcs(dframe, destname, sort_tied, 1)
        make_krona(destname)

    elif args.function == "assume_human":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        args = parser.parse_args(sub_args)
        dframe = pd.read_csv(args.input, sep="\t")
        filename = (args.input.split("/")[-1]).split(".")[0]
        destname = f"{output_dir}{filename}_output.txt"
        assume_human(dframe, destname)

    elif args.function == "count_taxid":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        args = parser.parse_args(sub_args)
        dframe = pd.read_csv(args.input, sep="\t")
        filename = (args.input.split("/")[-1]).split(".")[0]
        destname = f"{output_dir}{filename}_output.txt"
        count_taxids(dframe, destname)

main()
