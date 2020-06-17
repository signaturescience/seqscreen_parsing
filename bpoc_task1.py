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


# # go thru the revised file and count each bpoc
# def parse_bpocs(filename, destfilename, txtfilename):

# 	# initialize our counters
# 	# total_rows = 0
# 	# bpoc_rows = 0

# 	# adhesion_count = 0
# 	# secretion_count = 0
# 	# host_cell_death_count = 0
# 	# antibiotic_count = 0
# 	# invasion_count = 0
# 	# evasion_count = 0
# 	# cytotoxicity_count = 0
# 	# degrade_ecm_count = 0
# 	# disable_organ_count = 0

# 	# counts = [adhesion_count, secretion_count, host_cell_death_count,
# 	# 		antibiotic_count, invasion_count, evasion_count,
# 	# 		cytotoxicity_count, degrade_ecm_count, disable_organ_count]

# 	counts = defaultdict(int)

# 	bpoc_names = ["adhesion", "secretion", "host_cell_death",
# 			"antibiotic", "invasion", "evasion",
# 			"cytotoxicity", "degrade_ecm", "disable_organ"]

# 	# grab from pandas
# 	header_names = ["query", "taxid", "go", 
# 			"multi_taxids_confidence", "go_id_confidence",
# 			"adhesion", "secretion", "host_cell_death",
# 			"antibiotic", "invasion", "evasion",
# 			"cytotoxicity", "degrade_ecm", "disable_organ",
# 			"size", "organism", "gene_name", "uniprot", "uniprot evalue"]

# 	# initialize our dictionaries
# 	adhesion_dict = {}
# 	secretion_dict = {}
# 	host_cell_death_dict = {}
# 	antibiotic_dict = {}
# 	invasion_dict = {}
# 	evasion_dict = {}
# 	cytotoxicity_dict = {}
# 	degrade_ecm_dict = {}
# 	disable_organ_dict = {}

# 	dicts = [adhesion_dict, secretion_dict, host_cell_death_dict, 
# 			antibiotic_dict, invasion_dict, evasion_dict,
# 			cytotoxicity_dict, degrade_ecm_dict, disable_organ_dict]

# 	elems = ["taxid", "organism", "gene_name", "uniprot", "uniprot evalue"]

# 	for this_dict in dicts:
# 		for elem in elems:
# 			this_dict[elem] = defaultdict(int)

# 	# read and parse the file

# 	# df = pd.read_csv(filename, delimiter='\t')

# 	# # writer = csv.writer(open(destfilename,'wb'))

# 	# for row in len(df.index):
# 	# 	bpoc = False
# 	# 	# if row["query"]=="query":
# 	# 	# 	writer.writerow([row])
# 	# 	if row["adhesion"]!='-' and row["query"]!="query":
# 	# 		total_rows+=1
# 	# 		for i in range(9):
# 	# 			if (int(row[bpoc_names[i]]) > 0):
# 	# 				bpoc = True
# 	# 				counts[i]+=1
# 	# 				for elem in elems:
# 	# 					dicts[i][elem][row[elem]]+=1
# 	# 		if (bpoc):
# 	# 			bpoc_rows+=1
				


# 	# df_out.to_csv(destfilename, delimiter='\t')


# 	writer = csv.writer(open(destfilename,'wb'))

# 	with open(filename) as tsvfile:
# 		with open(destfilename, "w") as desttsvfile:
# 			reader = csv.DictReader(tsvfile, dialect='excel-tab')
# 			writer = csv.DictWriter(desttsvfile, fieldnames=header_names, delimiter='\t')
# 			writer.writeheader()
# 			for row in reader:
# 				bpoc = False
# 				if row["adhesion"]!='-' and row["query"]!="query":
# 					total_rows+=1
# 					for i in range(9):
# 						if (int(row[bpoc_names[i]]) > 0):
# 							bpoc = True
# 							counts[i]+=1
# 							for elem in elems:
# 								dicts[i][elem][row[elem]]+=1
# 					if (bpoc):
# 						bpoc_rows+=1
# 						writer.writerow(row)




# 	# write to the text file
# 	txt_file = open(txtfilename, "w")
# 	txt_file.write("Total rows: " + str(total_rows) + "\n")
# 	txt_file.write("Bpoc rows: " + str(bpoc_rows) + "\n")
# 	txt_file.write("Percentage of bpoc rows: " + str((1.0*bpoc_rows)/total_rows) + "\n")

# 	txt_file.write("\n\n")

# 	for i in range(len(counts)):
# 		txt_file.write("Percentage of " + bpoc_names[i] + ": " + str((1.0*counts[i])/total_rows) + "\n")

# 	txt_file.write("\n\n")
	
# 	for i in range(len(dicts)):
# 		for elem in elems:
# 			txt_file.write(bpoc_names[i] + " " + elem + " ")
# 			for j in dicts[i][elem]:
# 				txt_file.write("key: " + str(j) + ", value: " + str(dicts[i][elem][j]))
# 				txt_file.write("\n")
# 			txt_file.write("\n")
# 		txt_file.write("\n")

"""

filepath = "/Users/winnieli/Documents/summer2020microbes/"
def parse_bpocs_filepath(filepath, filename, destfilename, txtfilename):
	parse_bpocs(filepath+filename, filepath+destfilename, filepath+txtfilename)

parse_bpocs_filepath(filepath, "S01_trim25_fast_seqscreen_report.tsv", "S01_revised.tsv", "S01.txt")
parse_bpocs_filepath(filepath, 'S02_trim25_fast_seqscreen_report.tsv', 'S02_revised.tsv', 'S02.txt')
parse_bpocs_filepath(filepath, 'S03_trim25_fast_seqscreen_report.tsv', 'S03_revised.tsv', 'S03.txt')
parse_bpocs_filepath(filepath, 'S04_trim25_fast_seqscreen_report.tsv', 'S04_revised.tsv', 'S04.txt')
parse_bpocs_filepath(filepath, 'SRR10903401_combined_trim25_fast_seqscreen_report.tsv', 'SRR401_revised.tsv', 'SRR401.txt')
parse_bpocs_filepath(filepath, 'SRR10903402_combined_trim25_fast_seqscreen_report.tsv', 'SRR402_revised.tsv', 'SRR402.txt')
# parse_bpocs_filepath(filepath, 'SRR10971381_combined_trim25_seqscreen_report.tsv', 'SRR381_revised.tsv', 'SRR381.txt')

"""
