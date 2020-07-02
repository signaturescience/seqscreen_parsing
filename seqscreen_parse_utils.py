# -*- coding: utf-8 -*-
"""
seqscreen_parse_utils.py

All functions used for parsing scripts in this repo.
"""
import copy
import subprocess
import pandas as pd

pd.set_option('display.max_colwidth', 10000)

def count_taxids(dframe, destname):
    """

    Parameters
    ----------
    inputname : pandas dataframe representation of seqscreen output.
    destname : path to output file (str).

    Returns
    -------
    None.

    Writes  file containing taxid frequencies.

    """
    # dframe = pd.read_csv(inputname, sep="\t")
    tid = "taxid"

    #count frequency of each taxid. first element does regular count,
    #second element does relative count
    inp = [dframe[tid].value_counts(), dframe[tid].value_counts(normalize=True)]
    summary = pd.concat(inp, 1)
    summary.loc['total'] = summary.sum(numeric_only=True, axis=0)
    summary = summary.reset_index()
    summary.columns = [tid, "count", "percentage"]
    summary.to_csv(destname, sep="\t", index=False)

def assume_human(dframe, destname):
    """

    Parameters
    ----------
    inputname : path to input file (str).
    destname : path to output file (str).

    Returns
    -------
    None.

    Edits taxid column to assume each sequence that contains the human taxid is
    a human sample.

    """
    dframe.loc[dframe["multi_taxids_confidence"].str.contains("9606:"),
               "taxid"] = 9606
    dframe.to_csv(destname, sep="\t", index=False)

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
        pd.set_option('max_colwidth', 100000)
        cell = cell.split("confidence")[1]

    taxids = {int(k): float(v) for k, v in [s.split(':') for s in cell.split(",")]}
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
    """

    Parameters
    ----------
    cell : a string or pandas series representing a cell in the
    multi_taxids_confidence column.

    thresh : the max. number of taxids with the same confidence level.

    Returns
    -------
    list
        first element is the taxid with max. confidence once tied taxids
        are removed. second element is edited multi_taxids_confidence cell for
        that row.

    """
    #removes all taxids of equally high confidence from a cell
    if not isinstance(cell, str):
        cell = cell.to_string()
        cell = cell.split("confidence")[1]

    #converts string into dictionary
    taxids = {int(k): float(v) for k, v in [s.split(':') for s in cell.split(",")]}
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

def parse_funcs(dframe, destname, func, attr):
    """

    Parameters
    ----------
    inputname : pandas dataframe of input file.
    func : either sort_conf or sort_tied
    attr : attributes needed to run func.

    """
    mtid = "multi_taxids_confidence"

    df2 = dframe[mtid]
    df2 = df2.to_frame()
    df2 = df2.apply(func, 1, result_type='expand', args=[attr])
    df2 = df2.rename(columns={0:"taxid", 1:"multi_taxids_confidence"})

    df2 = pd.concat([dframe["query"], df2["taxid"], dframe["go"], df2["multi_taxids_confidence"],
                     dframe.iloc[:, 4:17]], 1)
    df2 = df2.dropna()
    df2["taxid"] = df2["taxid"].astype(int)
    df2.to_csv(destname, sep="\t", index=False)

def make_krona(infile):
    """

    Parameters
    ----------
    infile :  path to input file (str).

    Returns
    -------
    None.

    Generates a .tsv file as a Krona input.
    Generates Krona chart using that file.

    """
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
                    tsvname, "-o", f"{input_prefix}krona.html"], check=True)
