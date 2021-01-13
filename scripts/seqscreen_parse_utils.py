# -*- coding: utf-8 -*-
"""
seqscreen_parse_utils.py

All functions used for parsing scripts in this repo.
"""
import copy
import os
import subprocess
import pandas as pd
import numpy as np


#why is this here?
pd.set_option('display.max_colwidth', 1000000)

table_elements = ["query","taxid", "organism", "gene_name", "uniprot", "uniprot evalue"]
bpocs = ["adhesion", "secretion", "host_cell_death",
         "antibiotic", "invasion", "evasion",
         "cytotoxicity", "degrade_ecm", "disable_organ"]

def krona_from_slice(slice,file_root, name='GO_term'):
    """"
    takes in pandas dataframe slice,  and file prefix and writes out ktImportTaxonomy
    """
    input_data = slice.to_csv(header=False, sep="\t")
    subprocess.run(["ktImportTaxonomy", "-n", name, "-q", "1", "-t", "2", "-s", "3", "-", "-o", file_root + ".html"],
                   check=True, input=input_data.encode('utf_8'))

def krona_plot(inputfilename):
    """
    Takes in an input file (.tsv) and an output directory,
    Creates krona plots for that tsv
    """
    temp = pd.read_csv(inputfilename, sep="\t")
    dataframe = temp["multi_taxids_confidence"].str.split(",")
    dataframe = pd.concat([temp["query"], dataframe], 1)
    dataframe = dataframe.explode("multi_taxids_confidence")
    if dataframe.empty:
        return "No Entries"
    if dataframe is None:
        return 1
    final = pd.concat([dataframe["query"], dataframe["multi_taxids_confidence"]
                       .str.split(":", expand=True)], 1)

    outfile_name = f"{inputfilename}_krona.txt"
    outfile = open(outfile_name, "w")
    final.to_csv(outfile, sep="\t", header=False, index=False)
    outfile.close()

    krona_outfile_name = f"{outfile_name}.html"
#    krona_outfile = open(krona_outfile_name, "w")
    subprocess.Popen(
        f"ktImportTaxonomy -q 1 -t 2 -s 3 {outfile_name} -o {krona_outfile_name}",
        shell=True
        ).wait()
#    krona_outfile.close()
    return(0)

def remove_blanks_from_dataframe(dataframe,column, blank_pattern = '-'):
    """

    All rows that are non dash, aka have a GO term assigned
    filter based on column name in function
    if blank pattern is not '-', then must be specified in function call
    """
    idx = dataframe[dataframe[column] == blank_pattern].index
    df_valid = dataframe.drop(idx)
    valid_rows = len(df_valid.index)
    return(df_valid, valid_rows)

def bpoc_parse(dataframe, filename, output_dir):
    """
    Creates a revised file with only bpoc lines and
    a summary file with analysis of bpocs
    """
    elems = table_elements

    total_rows = len(dataframe.index)
    (df_valid, valid_rows) = remove_blanks_from_dataframe(dataframe, 'go')

    # All rows that have at least one BPoC
    df_bpocs = dataframe[dataframe[bpocs].replace('-', 0).astype(int).sum(1) > 0]
    bpoc_rows = len(df_bpocs.index)
    df_bpocs.to_csv(os.path.join(output_dir, filename + "_revised.tsv"), sep='\t', index=False)

    # What taxid, organism, gene_name, uniprot, and uniprot evalues
    # were assigned to the BPoCs within the sample?
    bpocs_series = df_bpocs[bpocs].astype(int).sum(0)
    bpoc_counts = pd.DataFrame(
        {'number':bpocs_series.values,
         'percentage':bpocs_series.values/valid_rows},
        index=bpocs_series.index)
    f_out = open(os.path.join(output_dir, filename + "_summary.txt"), "w") # summary file
    f_out.write(bpoc_counts.to_string())
    f_out.write("\n")
    f_out.write(f"Total rows including dashes: {total_rows}")
    f_out.write("\n")
    f_out.write(f"Total rows excluding dashes: {valid_rows}")
    f_out.write("\n")
    f_out.write("Percentage: bpoc in sample/number of rows without dashes: ")
    f_out.write(str(bpoc_rows/(1.0*valid_rows)))
    f_out.write("\n\n")

    # What taxid, organism, gene_name, uniprot, and uniprot
    # evalues were assigned to the BPoCs within the sample?
    f_out.write("Taxid, organism, gene_name, uniprot, and uniprot evalues per BPoC:")
    for bpoc in bpocs:
        if bpoc_counts.at[bpoc, 'number'] > 0:
            bpoc_elems = df_bpocs[df_bpocs[bpoc].astype(int) > 0][elems]
            f_out.write(f"\n{bpoc} \n")
            elems_in_bpoc = pd.Series(bpoc_elems.transpose().to_numpy().tolist(),
                                      index=bpoc_elems.columns).apply(lambda x: ', '.join(set(x)))
            #ideally this would not be a for loop,
            #but dislike formatting for elems_in_bpoc.to_string()
            records = elems_in_bpoc.csv()
            f_out.write(records)
            #for index, value in elems_in_bpoc.items():
            #    f_out.write(f"{index}: {value} \n")

    f_out.close()

## Takes in godag, dataframe and list of GO terms, returns dataframe of only query terms associated with GO terms in list.
def parse_GO_terms(godag, dataframe, go_nums):
    return_df = pd.DataFrame(columns=['GO_term', 'query', 'organism', 'associated_GO_terms', 'multi_taxids_confidence',
                                      'taxid', 'gene_name', 'uniprot', 'uniprot evalue'], dtype='str')
    iter = 0
    """
    This structure is horribly inefficient.  Needs to be refactored to do the following:
    generate expanded dataframe of ALL queries and ALL GO terms
    Filter to only include those rows that contain GO term
    check query/GO combo for any GO term children, place in associated_GO_terms
    """
    rows_raw = len(dataframe.index)
    dataframe2 = dataframe
    (dataframe, num_rows) = remove_blanks_from_dataframe(dataframe, 'go')
    print(num_rows, "of", rows_raw, "with annotated GO terms kept")
    dataframe['go_id_confidence'] = dataframe['go_id_confidence'].str.split(";")
    temp_taxID = dataframe['multi_taxids_confidence'].str.split(';')
    expanded_dataframe = dataframe.explode('go_id_confidence')
    expanded_dataframe['go'] = expanded_dataframe['go_id_confidence'].str.replace("\[.*", "", regex=True)
    godag_keys = godag.keys()
    total_iter = 0
    for go in go_nums:
        print(go)
        string_iter = 0
        slice = expanded_dataframe.loc[expanded_dataframe['go'] == go]
        go_family = [go]
        if go in godag_keys:
            for value in list(godag[go].get_all_children()):
                go_family.append(value)
        # slice2 = expanded_dataframe.loc[expanded_dataframe['go'].isin(go_family)]
        if len(slice.index) > 0:
            print("0", "/", len(slice.index))
            queries = list(set(slice['query']))
            string_iter = 0
            for query in queries:
                total_iter += 1
                string_iter += 1
                if string_iter % 1000 == 0:
                    print(string_iter, "/", len(queries))
                slice2 = dataframe2.loc[dataframe2['query'] == query]
#                query_list = expanded_dataframe.loc[expanded_dataframe['query'] == query, 'go_id_confidence']
#                query_list = np.intersect1d(query_list, go_family)
                fl = expanded_dataframe['go']['query' == query and 'go' in go_family]
                query_list=np.intersect1d(fl,go_family)
                idx = slice2.index[0]
                return_df.loc[total_iter] = {'GO_term': go, 'query': query,
                                             'organism': slice2['organism'][idx],
                                             'associated_GO_terms': ";".join(query_list), 'multi_taxids_confidence': slice2['multi_taxids_confidence'], 'taxid': slice2['taxid'][idx],
                                             'gene_name': slice2['gene_name'][idx], 'uniprot': slice2['uniprot'][idx],
                                             'uniprot evalue': slice2['uniprot evalue'][idx]}  # dicer
        print(string_iter, "/", string_iter, ":", go, "Complete")
    return_df.to_csv("test.csv")
    return (return_df)


def collapse_GO_results(dataframe, group_by, count):
    return dataframe.groupby(group_by)[count].count()



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
        cell = cell.split("confidence")[1]

    #turns cell into dictionary
    taxids = {int(k): float(v) for k, v in [s.split(':') for s in cell.split(",")]}
    final_tids = copy.deepcopy(taxids)

    #removes taxids below conf
    for key, val in taxids.items():
        if val < conf:
            final_tids.pop(key)

    if final_tids != {}:
        strdict = repr(taxids).replace("{", "").replace("}", "")
        return [max(final_tids, key=final_tids.get), strdict]
    return [None, None]

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
        strdict = repr(taxids).replace("{", "").replace("}", "")
        return [max(taxids, key=taxids.get), strdict]
    return [None, None]

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
    df2 = df2.rename({0:"taxid", 1:"multi_taxids_confidence"}, axis='columns')
    df2 = pd.concat([dframe["query"], df2["taxid"], dframe["go"], df2["multi_taxids_confidence"],
                     dframe.iloc[:, 4:17]], 1)
    df2 = df2.dropna()
    df2["taxid"] = df2["taxid"].astype(int)
    if len(df2.index) < 1:
        print("no samples in input file fit given criteria, empty file returned")
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

    if len(temp.index) < 1:
        return

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


def create_output_directory(directory_name):
    if not os.path.exists(directory_name):
        try:
            os.mkdir(directory_name)
        except:
            print("could not create directory", directory_name)
            return 1
    return 0
