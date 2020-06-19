from collections import defaultdict
import csv
import pandas as pd
import argparse
import pathlib

# parse the inputted .tsv file (needs to be full path)
parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str, help="input a .tsv file")
args = parser.parse_args()
filename = pathlib.Path(args.input_file)


# remove rows where there are no bpocs, create a revised file
df = pd.read_csv(filename, sep='\t', dtype=str)

# filename_short = str(filename).split('/')[-1].split('.')[-2]

bpocs = ["adhesion", "secretion", "host_cell_death",
			"antibiotic", "invasion", "evasion",
			"cytotoxicity", "degrade_ecm", "disable_organ"]
elems = ["taxid","organism","gene_name","uniprot", "uniprot evalue"]

df_bpocs = df[df[bpocs].replace('-',0).astype(int).sum(1)>0]
df_bpocs.to_csv(f"{filename}_revised.tsv", sep='\t', index = False)

# What taxid, organism, gene_name, uniprot, and uniprot evalues 
# were assigned to the BPoCs within the sample?
f_out = open(f"{filename}_summary.txt", "w")
bpoc_counts = df_bpocs[bpocs].astype(int).sum(0)
f_out.write(bpoc_counts.to_string())
f_out.write("\n")
f_out.write(f"Percentage of bpoc in sample:{len(df_bpocs.index)/len(df.index)}")
f_out.write("\n\n")

# TODO: process multi taxids
# What taxid, organism, gene_name, uniprot, and uniprot evalues were assigned to the BPoCs within the sample?
f_out.write("Taxid, organism, gene_name, uniprot, and uniprot evalues per BPoC:")
for bpoc in bpocs:
    if(bpoc_counts[bpoc]>0):
        bpoc_elems = df_bpocs[df_bpocs[bpoc].astype(int)>0][elems]
        f_out.write(f"\n{bpoc} \n")
        elems_in_bpoc = pd.Series(bpoc_elems.transpose().to_numpy().tolist(),
                                  index=bpoc_elems.columns).apply(lambda x: ', '.join(set(x)))
        f_out.write(f"{elems_in_bpoc.to_string()} \n")
f_out.close()






