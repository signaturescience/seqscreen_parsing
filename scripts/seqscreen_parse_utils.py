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

def bpoc_parse(dataframe, filename, output_dir):
    """
    Creates a revised file with only bpoc lines and
    a summary file with analysis of bpocs
    """
    bpocs = ["adhesion", "secretion", "host_cell_death",
             "antibiotic", "invasion", "evasion",
             "cytotoxicity", "degrade_ecm", "disable_organ"]

    elems = ["taxid", "organism", "gene_name", "uniprot", "uniprot evalue"]

    total_rows = len(dataframe.index)

    # All rows that are non dash, aka have a GO term assigned
    idx = dataframe[dataframe['go'] == '-'].index
    df_valid = dataframe.drop(idx)
    valid_rows = len(df_valid.index)

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

## Takes in godag, dataframe and list of GO terms, returns dataframe of only query terms associated with GO terms in list.
def parse_GO_terms(godag,dataframe, go_nums):
    return_df = pd.DataFrame(columns=['GO_term','query','organism','associated_GO_terms','taxid', 'gene_name', 'uniprot', 'uniprot evalue'])
    iter=0
    for go in go_nums:
        #    family_tree = [go, godag[go].get_all_children()]
        family_tree = [go]
        blah = []
        for j in family_tree:
            if type(j) == str:
                blah.append(j)
            else:
                blah.extend(j)
        family_tree = blah
        for i, row in dataframe.iterrows():
            input_go_terms = row['go'].split(";")
            input_confidence = row['go_id_confidence'].split(";")
            list = []
            for j in range(0, (len(input_go_terms) - 1)):
                term = input_go_terms[j]
                confidence = input_confidence[j]
                if term == go:
                    # print('go term', go, 'found in', i, "as", input_go_terms)
                    list.append(confidence)
                #                break
                elif term in godag.keys():
                    if go in godag[term].get_all_parents():
                        list.append(confidence)
            if len(list) > 0:
                iter += 1
                return_df.loc[str(iter)] = pd.Series({'GO_term':go,'query':row['query'],'organism':row['organism'],'associated_GO_terms':list,'taxid':row['taxid'], 'gene_name':row['gene_name'],'uniprot':row['uniprot'],'uniprot evalue':row['uniprot evalue']})
                #print(go,row['query'])
    return(return_df)

def collapse_GO_results(dataframe):
    return dataframe.groupby(['GO_term','taxid'])['taxid'].count()

'''
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
'''

pd.set_option('display.max_colwidth', 1000000)

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


def output_directory(directory_name):
    if not os.path.exists(directory_name):
        try:
            os.mkdir(directory_name)
        except:
            print("could not create directory", directory_name)
            return 1
    return 0
