"""
REU_2020 microbes task 2
https://trello.com/b/tDbNjSfc/reu2020

@author: wwl3
"""
import argparse
import os
import pathlib
import pandas as pd

def get_tied_taxids(multi_taxids_confidence):
    """
    For each multi_taxid_confidence cell, gets a list of the top taxids
    that have tied confidence levels
    """
    taxids = []
    data = multi_taxids_confidence.split(",")
    maximum = 0
    for dat in data:
        taxid_conf = dat.split(":")
        if float(taxid_conf[1]) > maximum:
            maximum = float(taxid_conf[1])
    for dat in data:
        taxid_conf = dat.split(":")
        if float(taxid_conf[1]) >= maximum:
            taxids.append(taxid_conf[0])
    return taxids

def go_term_parse(dataframe, go_num, filename, output_dir):
    """
    Creates a revised file with only lines with the go term specified
    and a summary file with analysis of that go term
    """

    columns = ["multi_taxids_confidence", "organism", "gene_name", "uniprot"]

    # create a new .tsv file with just rows with that go number
    df_go = dataframe[dataframe["go"].apply(lambda x: go_num in x)]
    df_go.to_csv(f"{output_dir}{filename}_go_num_{go_num}_revised.tsv", sep='\t', index=False)

    # create a text file  of summary statistics for that go number
    f_out = open(f"{output_dir}{filename}_{go_num}_summary.txt", "w") # summary file
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

def main():
    """
    Runs the go term parsing from the command line
    takes in a .tsv file and a list of go terms as input
    """

    # parse the inputted .tsv file
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="input a .tsv file")
    parser.add_argument(
        "go_numbers", type=str, help="input go numbers separated by commas, with no spaces")
    args = parser.parse_args()
    filename = pathlib.PurePath(args.input_file).stem
    output_dir = "outputs/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # remove rows where there are no bpocs, create a revised file
    dataframe = pd.read_csv(pathlib.Path(args.input_file), sep='\t', dtype=str)
    go_nums = args.go_numbers.split(',')
    for go_num in go_nums:
        go_term_parse(dataframe, go_num, filename, output_dir)

main()
