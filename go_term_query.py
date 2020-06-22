"""
REU_2020 microbes task 2
https://trello.com/b/tDbNjSfc/reu2020

@author: wwl3
"""

# What visualizations could be created?
import pathlib
import argparse
import pandas as pd

# parse the inputted .tsv file (needs to be full path)
parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="input a .tsv file")
parser.add_argument(
    "go_numbers", type=str, help="input go numbers separated by semicolons, with no spaces")
args = parser.parse_args()
filename = pathlib.Path(args.input_file)
go_nums = args.go_numbers.split(';')


# remove rows where there are no bpocs, create a revised file
df = pd.read_csv(filename, sep='\t', dtype=str)

# filename_short = str(filename).split('/')[-1].split('.')[-2]
columns = ["taxid", "organism", "gene_name", "uniprot"]

for go_num in go_nums:
    # create a new .tsv file with just rows with that go number
    df_go = df[df["go"].apply(lambda x: go_num in x)]
    df_go.to_csv(f"{filename}_go_num_{go_num}_revised.tsv", sep='\t', index=False)

    # create a text file  of summary statistics for that go number
    f_out = open(f"{filename}_go_num_{go_num}_summary.txt", "w")
    f_out.write(f"Information associated with go term {go_num} for sample: {filename}")
    f_out.write("\n\n")
    f_out.write(
        f"Percentage of rows with go term {go_num} in sample:{len(df_go.index)/len(df.index)}")
    f_out.write("\n\n")

    # summarize the proteins and organisms
    for col in columns:
        # if col == "taxid":
        #   df_go_counts = df_go.groupby("multi_taxids_confidence").
        #   .apply({
        #       lambda x:
  #               taxids = []
        #       data = x.split(",")
        #       maximum = 0
        #       for dat in data:
        #           taxid_conf = dat.split(":")
        #           if taxid_conf[1] > maximum:
        #               maximum = taxid_conf[1]
        #       for dat in data:
        #           taxid_conf = dat.split(":")
        #           if taxid_conf[1] >= maximum:
        #               taxixds.append(taxid_conf[0]))
        #   .count().loc[:, "query"].sort_values(ascending=False)
            # count().loc[:, "query"].sort_values(ascending=False)
        df_go_counts = df_go.groupby(col).count().loc[:, "query"].sort_values(ascending=False)
        f_out.write(f"{df_go_counts.to_string()} \n\n")
    f_out.close()
