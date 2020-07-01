# -*- coding: utf-8 -*-
"""
seqscreen_parse_utils.py

All functions used for parsing scripts in this repo.
"""
import copy
import os
import subprocess
import pandas as pd

def krona_plot(inputfilename):
    """
    Takes in an input file (.tsv) and an output directory,
    Creates krona plots for that tsv
    """
    temp = pd.read_csv(inputfilename, sep="\t")
    dataframe = temp["multi_taxids_confidence"].str.split(",")
    dataframe = pd.concat([temp["query"], dataframe], 1)
    dataframe = dataframe.explode("multi_taxids_confidence")
    final = pd.concat([dataframe["query"], dataframe["multi_taxids_confidence"]
                       .str.split(":", expand=True)], 1)

    outfile_name = f"{inputfilename}_krona.txt"
    outfile = open(outfile_name, "w")
    final.to_csv(outfile, sep="\t", header=False, index=False)

    krona_outfile_name = f"{outfile_name}.html"
    open(krona_outfile_name, "w")
    subprocess.Popen(f"ktImportTaxonomy -q 1 -t 2 -s 3 {outfile_name} -o {krona_outfile_name}", shell=True).wait()

def bpoc_parse(dataframe, filename, output_dir):
    """
    Creates a revised file with only bpoc lines and
    a summary file with analysis of bpocs
    """
    bpocs = ["adhesion", "secretion", "host_cell_death",
             "antibiotic", "invasion", "evasion",
             "cytotoxicity", "degrade_ecm", "disable_organ"]

    elems = ["taxid", "organism", "gene_name", "uniprot", "uniprot evalue"]

    df_bpocs = dataframe[dataframe[bpocs].replace('-', 0).astype(int).sum(1) > 0]
    df_bpocs.to_csv(os.path.join(output_dir, filename + "_revised.tsv"), sep='\t', index=False)

    # What taxid, organism, gene_name, uniprot, and uniprot evalues
    # were assigned to the BPoCs within the sample?
    bpoc_counts = df_bpocs[bpocs].astype(int).sum(0)
    f_out = open(os.path.join(output_dir, filename + "_summary.txt"), "w") # summary file
    f_out.write(bpoc_counts.to_string())
    f_out.write("\n")
    f_out.write(f"Percentage of bpoc in sample:{len(df_bpocs.index)/len(dataframe.index)}")
    f_out.write("\n\n")

    # What taxid, organism, gene_name, uniprot, and uniprot
    # evalues were assigned to the BPoCs within the sample?
    f_out.write("Taxid, organism, gene_name, uniprot, and uniprot evalues per BPoC:")
    for bpoc in bpocs:
        if bpoc_counts[bpoc] > 0:
            bpoc_elems = df_bpocs[df_bpocs[bpoc].astype(int) > 0][elems]
            f_out.write(f"\n{bpoc} \n")
            elems_in_bpoc = pd.Series(bpoc_elems.transpose().to_numpy().tolist(),
                                      index=bpoc_elems.columns).apply(lambda x: ', '.join(set(x)))
            #ideally this would not be a for loop,
            #but dislike formatting for elems_in_bpoc.to_string()
            for index, value in elems_in_bpoc.items():
                f_out.write(f"{index}: {value} \n")

    f_out.close()

def get_tied_taxids(multi_taxids_confidence):
    """
    For each multi_taxid_confidence cell, gets a list of the top taxids
    that have tied confidence levels
    """
    taxids = []
    data = multi_taxids_confidence.split(",")
    maximum = max(float(cell.split(':')[1]) for cell in data)
    taxids = [cell.split(':')[0] for cell in data if float(cell.split(':')[1]) == maximum]

    return taxids

def go_term_parse(dataframe, go_num, filename, output_dir):
    """
    Creates a revised file with only lines with the go term specified
    and a summary file with analysis of that go term
    """

    columns = ["multi_taxids_confidence", "organism", "gene_name", "uniprot"]

    # create a new .tsv file with just rows with that go number
    df_go = dataframe[dataframe["go"].apply(lambda x: go_num in x)]
    df_go.to_csv(os.path.join(output_dir, filename + f"_go_num_{go_num}_revised.tsv"),
                 sep='\t', index=False)

    # create a text file  of summary statistics for that go number
    f_out = open(os.path.join(output_dir, filename + f"_{go_num}_summary.txt"), "w")
    f_out.write(f"Information associated with go term {go_num} for sample: {filename}")
    f_out.write("\n\n")
    f_out.write(f"Percentage of rows with go term {go_num} in ")
    f_out.write(f"sample: {len(df_go.index)/len(dataframe.index)}")
    f_out.write("\n\n")

    # summarize the proteins and organisms
    for col in columns:
        if col == "multi_taxids_confidence":
            df_go_counts = df_go[col].apply(get_tied_taxids).explode().value_counts()
        else:
            df_go_counts = df_go.groupby(col).count().loc[:, "query"].sort_values(ascending=False)
        f_out.write(f"{col} \n")
        f_out.write(f"{df_go_counts.to_string()} \n\n")
    f_out.close()
    return f_out

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
