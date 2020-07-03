#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
taxid_parse.py

Script to extract information regarding taxon abundances in sequences.
"""
import argparse
import os
import pandas as pd
import seqscreen_parse_utils as spu

def main():
    """
    Runs taxid parsing from command line.
    Takes in .tsv file as input.

    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("input", type=str, help="filepath to input .tsv")
    parser.add_argument("--parse_conf", type=float,
                        help="remove taxids below input confidence")
    parser.add_argument("--thresh_tied", type=int,
                        help="remove tied taxids above given threshold")
    parser.add_argument("--all_tied",
                        help="remove all tied taxids")
    parser.add_argument("--assume_human",
                        help="assume human sequence if human taxid is in multi_taxid")
    parser.add_argument("--count_taxid",
                        help="count frequency and percentage of multi_taxids")

    args = parser.parse_args()
    dframe = pd.read_csv(args.input, sep="\t")

    tempfile = args.input.split(".")[0].split("/")[-1]
    dirname = args.input.split(tempfile)[0] + "outputs/"

    if not os.path.exists(dirname):
        os.makedirs(dirname)
    destname = dirname + tempfile + ".tsv"

    if args.all_tied:
        spu.parse_funcs(dframe, destname, spu.sort_tied, 1)
        spu.make_krona(destname)

    if args.parse_conf is not None:
        spu.parse_funcs(dframe, destname, spu.sort_conf, args.parse_conf)
        spu.make_krona(destname)

    if args.thresh_tied is not None:
        spu.parse_funcs(dframe, destname, spu.sort_tied, args.thresh_tied)
        spu.make_krona(destname)

    if args.assume_human:
        spu.assume_human(dframe, destname)

    if args.count_taxid:
        spu.count_taxids(dframe, destname)

main()
