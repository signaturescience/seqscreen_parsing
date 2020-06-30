#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:28:04 2020

@author: Student
"""
import argparse
import copy
import subprocess
import pandas as pd

def count_taxids(inputname):
    dframe = pd.read_csv(inputname, sep="\t")
    tid = "taxid"

    #count frequency of each taxid. first element does regular count,
    #second element does relative count
    inp = [dframe[tid].value_counts(), dframe[tid].value_counts(normalize=True)]
    summary = pd.concat(inp, 1)
    summary.loc['total'] = summary.sum(numeric_only=True, axis=0)
    summary = summary.reset_index()
    summary.columns = [tid, "count", "percentage"]
    name = inputname.split(".")[0]
    summary.to_csv(f"{name}output.txt", sep="\t", index=False)

def assume_human(inputname):
    dframe = pd.read_csv(inputname, sep="\t")
    dframe.loc[dframe["multi_taxids_confidence"].str.contains("9606:"),
               "taxid"] = 9606
    name = inputname.split(".")[0]
    dframe.to_csv(f"{name}output.txt", sep="\t", index=False)

def sort_conf(cell, conf):
    """

    Parameters
    ----------
    cell : string or pandas series representing contents of cell
    conf : float representing confidence level

    Returns
    -------
    list
        taxid with max confidence level, updated string of taxid confidences

    """
    if not isinstance(cell, str):
        cell = cell.to_string()
        cell = cell.split("confidence")[1]

    taxids = eval("{" + cell + "}")
    final_tids = copy.deepcopy(taxids)

    #removes taxids below conf
    for key, val in taxids.items():
        if val < conf:
            final_tids.pop(key)

    if final_tids != {}:
        temp = repr(final_tids).replace("{", "")
        strdict = temp.replace("}", "")
        return [max(final_tids, key=final_tids.get), strdict]

def sort_tied(cell, thresh):
    #removes all taxids of equally high confidence from a cell
    if not isinstance(cell, str):
        cell = cell.to_string()
        cell = cell.split("confidence")[1]

    #converts string into dictionary
    taxids = eval("{" + cell + "}")
    rev = {}

    for key, value in taxids.items():
        rev.setdefault(value, set()).add(key)

    for key in rev:
        if len(rev[key]) > thresh:
            for tid in rev[key]:
                taxids.pop(tid)

    if taxids != {}:
        temp = repr(taxids).replace("{", "")
        strdict = temp.replace("}", "")
        return [max(taxids, key=taxids.get), strdict]

def parse_funcs(inputname, func, attr):
    """

    Parameters
    ----------
    inputname : path to input file
    func : either sort_conf or sort_tied
    attr : attributes needed to run func.

    """
    dframe = pd.read_csv(inputname, sep="\t")
    mtid = "multi_taxids_confidence"

    df2 = dframe[mtid]
    df2 = df2.to_frame()
    df2 = df2.apply(func, 1, result_type='expand', args=[attr])
    df2 = df2.rename(columns={0:"taxid", 1:mtid})

    df2 = pd.concat([dframe["query"], df2["taxid"], dframe["go"], df2[mtid],
                     dframe.iloc[:, 4:17]], 1)
    df2 = df2.dropna()
    df2["taxid"] = df2["taxid"].astype(int)

    name = inputname.split(".")[0]
    df2.to_csv(f"{name}output.txt", sep="\t", index=False)
    return f"{name}output.txt"

def make_krona(infile):
    temp = pd.read_csv(infile, sep="\t")
    mtid = "multi_taxids_confidence"
    dframe = temp[mtid].str.split(",")
    dframe = pd.concat([temp["query"], dframe], 1)
    dframe = dframe.explode(mtid)
    final = pd.concat([dframe["query"],
                       dframe[mtid].str.split(":", expand=True)], 1)
    input_prefix = infile.split(".")[0]
    tsvname = f"{input_prefix}_kt_in.tsv"
    final.to_csv(tsvname, sep="\t", header=False, index=False)
    subprocess.run(["ktImportTaxonomy", "-q", "1", "-t", "2", "-s", "3",
                    tsvname, "-o", f"{input_prefix}krona.html"])

def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("function",
                        nargs="?",
                        choices=['parse_conf', 'thresh_tied', 'all_tied',
                                 'assume_human', "count_taxid"],)

    args, sub_args = parser.parse_known_args()

    if args.function == "parse_conf":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        parser.add_argument('-c', '--confidence', type=float)
        args = parser.parse_args(sub_args)
        fname = parse_funcs(args.input, sort_conf, args.confidence)
        make_krona(fname)

    elif args.function == "thresh_tied":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        parser.add_argument('-t', '--threshold', type=int)
        args = parser.parse_args(sub_args)
        fname = parse_funcs(args.input, sort_tied, args.threshold)
        make_krona(fname)

    elif args.function == "all_tied":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        args = parser.parse_args(sub_args)
        fname = parse_funcs(args.input, sort_tied, 1)
        make_krona(fname)

    elif args.function == "assume_human":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        args = parser.parse_args(sub_args)
        assume_human(args.input)

    elif args.function == "count_taxid":
        parser = argparse.ArgumentParser()
        parser.add_argument("input", type=str)
        args = parser.parse_args(sub_args)
        count_taxids(args.input)

main()
    