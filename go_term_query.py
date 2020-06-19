# End users interested in specific GO terms would like to retrieve 
#a list of all the rows that contain a specific GO term ID, or a 
#list of GO term IDs of interest to them. It would be useful to 
#summarize the proteins and organisms associated with each of the 
#GO terms of interest in a SeqScreen output file.

# How could a capability be setup for users to query specific 
# GO terms within a sample, and return all of the information associated 
#with each of those GO terms (e.g., organism, gene_name, uniprot)?

# What visualizations could be created?
from collections import defaultdict
import csv
import pandas as pd
import argparse
import pathlib

# parse the inputted .tsv file (needs to be full path)
parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="input a .tsv file")
parser.add_argument("go_numbers", type=str, help="input go numbers separated by semicolons, with no spaces")
args = parser.parse_args()
filename = pathlib.Path(args.input_file)
go_nums = args.go_numbers.split(';')


# remove rows where there are no bpocs, create a revised file
df = pd.read_csv(filename, sep='\t', dtype=str)

# filename_short = str(filename).split('/')[-1].split('.')[-2]
columns = ["taxid", "organism", "gene_name"]

for go_num in go_nums:
	# create a new .tsv file with just rows with that go number
	df_go = df[df["go"].apply(lambda x: go_num in x)]
	df_go.to_csv(f"{filename}_go_num_{go_num}_revised.tsv", sep='\t', index = False)

	# create a text file  of summary statistics for that go number
	f_out = open(f"{filename}_go_num_{go_num}_summary.txt", "w")
	f_out.write(f"Information associated with go term {go_num} for sample: {filename}")
	f_out.write("\n\n")
	f_out.write(f"Percentage of rows with go term {go_num} in sample:{len(df_go.index)/len(df.index)}")
	f_out.write("\n\n")
	# summarize the proteins and organisms
	for col in columns:
		df_go_counts = df_go.groupby(col).count().loc[:, "query"].sort_values(ascending=False)
		f_out.write(f"{df_go_counts.to_string()} \n\n")
	f_out.close()

