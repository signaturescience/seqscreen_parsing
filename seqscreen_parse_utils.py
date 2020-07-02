# -*- coding: utf-8 -*-
"""
seqscreen_parse_utils.py

All functions used for parsing scripts in this repo.
"""
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
    outfile.close()

    krona_outfile_name = f"{outfile_name}.html"
    krona_outfile = open(krona_outfile_name, "w")
    subprocess.Popen(f"ktImportTaxonomy -q 1 -t 2 -s 3 {outfile_name} -o {krona_outfile_name}", shell=True).wait()
    krona_outfile.close()

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
