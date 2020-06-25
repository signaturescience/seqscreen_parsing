#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:28:04 2020

@author: Student
"""
import copy
import pandas as pd
from argparse import ArgumentParser

def krona_input(infile, outfile): 
    temp = pd.read_csv(infile, sep="\t")
    df = temp["multi_taxids_confidence"].str.split(",")
    df = pd.concat([temp["query"], df], 1)
    df = df.explode("multi_taxids_confidence")
    final = pd.concat([df["query"], df["multi_taxids_confidence"].str.split(":", expand=True)], 1)
    final.to_csv(outfile, sep="\t", header=False, index=False)

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

    if final_tids == {}:
        return None
    else:
        temp = repr(final_tids).replace("{", "")
        strdict = temp.replace("}", "")
        return [max(final_tids, key=final_tids.get), strdict]
        # return str(max(final_tids, key=final_tids.get)) + "; " + strdict

def parse_funcs(inputname, func, attr):
    """

    Parameters
    ----------
    inputname : path to input file
    func : either sort_conf or sort_tied
    attr : attributes needed to run func.

    """
    df = pd.read_csv(inputname, sep="\t")
    l = "multi_taxids_confidence"

    b = df[l]
    b = b.to_frame()
    b = b.apply(func, 1, result_type='expand', args=[attr])
    b = b.rename(columns={0:"taxid", 1:l})

    b = pd.concat([df["query"], b["taxid"], df["go"], b[l], df.iloc[:, 4:17]], 1)
    b = b.dropna()
    b["taxid"] = b["taxid"].astype(int)

    n = inputname.split(".")[0]
    b.to_csv(f"{n}output.txt", sep="\t", index=False)

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

    for key in rev.keys():
        if len(rev[key]) > thresh:
            for tid in rev[key]:
                taxids.pop(tid)
    if taxids == {}:
        return None
    else:
        temp = repr(taxids).replace("{", "")
        strdict = temp.replace("}", "")
        return [max(taxids, key=taxids.get), strdict]

def count_taxids(inputname, destname):
    df = pd.read_csv(inputname, sep="\t")
    l = "taxid"

    #count frequency of each taxid. first element does regular count,
    #second element does relative count
    inp = [df[l].value_counts(), df[l].value_counts(normalize=True)]
    summary = pd.concat(inp, 1)
    summary.loc['total'] = summary.sum(numeric_only=True, axis=0)
    summary = summary.reset_index()
    summary.columns = ["taxid", "count", "percentage"]
    summary.to_csv(destname, sep="\t", index=False)

def assume_human(inputname):
    df = pd.read_csv(inputname, sep="\t")
    l = "multi_taxids_confidence"
    df.loc[df[l].str.contains("9606:"), "taxid"] = 9606
    n = inputname.split(".")[0]
    df.to_csv(f"{n}output.txt", sep="\t", index=False)
    
def parse_conf(inputname, attrs):
    parse_funcs(inputname, sort_conf, attrs)
    return "done"

def remove_tied(inputname, attrs):
    parse_funcs(inputname, sort_tied, attrs)
    return "done"

def remove_equal(inputname):
    parse_funcs(inputname, sort_tied, 1)
    return "done"



def main():
    parser = ArgumentParser(description="calls desired function")
    parser.add_argument("function",
                        nargs="?",
                        choices=['parse_conf', 'remove_tied', 'remove_equal',
                                 'assume_human', "count_taxids"],
                        default='parse_conf',)
    args, sub_args = parser.parse_known_args()

    if args.function == "parse_conf":
        parser = ArgumentParser()
        parser.add_argument('-i', '--inputname')
        parser.add_argument('-d', '--destname')
        parser.add_argument('-c', '--confidence')
        args = parser.parse_args(sub_args)
        parse_conf(args.inputname, args.confidence)

    if args.function == "remove_tied":
        parser = ArgumentParser()
        parser.add_argument('-i', '--inputname')
        parser.add_argument('-d', '--destname')
        parser.add_argument('-t', '--threshold')
        args = parser.parse_args(sub_args)
        remove_tied(args.inputname, args.threshold)

    elif args.function == ("remove_equal" or "assume_human" or "count_taxids"):
        parser = ArgumentParser()
        parser.add_argument('-i', '--inputname')
        parser.add_argument('-d', '--destname')
        args = parser.parse_args(sub_args)
        args.function(args.inputname, args.destname)

main()
