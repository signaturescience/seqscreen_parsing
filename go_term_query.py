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

def go_query(go_numbers, filename, destfilename, txtfilename):
	# initialize our counters
	total_rows = 0
	go_rows = 0

	adhesion_count = 0
	secretion_count = 0
	host_cell_death_count = 0
	antibiotic_count = 0
	invasion_count = 0
	evasion_count = 0
	cytotoxicity_count = 0
	degrade_ecm_count = 0
	disable_organ_count = 0

	counts = [adhesion_count, secretion_count, host_cell_death_count,
			antibiotic_count, invasion_count, evasion_count,
			cytotoxicity_count, degrade_ecm_count, disable_organ_count]

	bpoc_names = ["adhesion", "secretion", "host_cell_death",
			"antibiotic", "invasion", "evasion",
			"cytotoxicity", "degrade_ecm", "disable_organ"]

	header_names = ["query", "taxid", "go", 
			"multi_taxids_confidence", "go_id_confidence",
			"adhesion", "secretion", "host_cell_death",
			"antibiotic", "invasion", "evasion",
			"cytotoxicity", "degrade_ecm", "disable_organ",
			"size", "organism", "gene_name", "uniprot", "uniprot evalue"]

	# initialize our dictionaries
	# taxid_dict = {}
	adhesion_dict = {}
	secretion_dict = {}
	host_cell_death_dict = {}
	antibiotic_dict = {}
	invasion_dict = {}
	evasion_dict = {}
	cytotoxicity_dict = {}
	degrade_ecm_dict = {}
	disable_organ_dict = {}
	organism_dict = {}
	# gene_name_dict = {}
	# uniprot_dict = {}
	# uniprot_evalue_dict = {}

	bpoc_dicts = [adhesion_dict, secretion_dict, host_cell_death_dict, 
			antibiotic_dict, invasion_dict, evasion_dict,
			cytotoxicity_dict, degrade_ecm_dict, disable_organ_dict]

	#elem_dicts = [taxid_dict, organism_dict, gene_name_dict, uniprot_dict, uniprot_evalue_dict]

	elems = ["taxid", "organism", "gene_name", "uniprot", "uniprot evalue"]

	for this_dict in bpoc_dicts:
		for elem in elems:
			this_dict[elem] = defaultdict(int)

	go_dict = {}
	for num in go_numbers:
		go_dict[str(num)] = {}
		for elem in elems:
			go_dict[str(num)][elem] = defaultdict(int)

	# read and parse the file
	with open(filename) as tsvfile:
		with open(destfilename, "w") as desttsvfile:
			reader = csv.DictReader(tsvfile, dialect='excel-tab')
			writer = csv.DictWriter(desttsvfile, fieldnames=header_names, delimiter = '\t')
			writer.writeheader()

			for row in reader:
				go_flags = set([])

				# if the row is valid (aka not header or row with dashes)
				if row["adhesion"]!='-' and row["query"]!="query":
					total_rows+=1
					this_golist = row["go"].split(";")
					for go_num1 in this_golist:
						for go_num2 in go_numbers:
							if go_num1 == go_num2:
								go_flags.add(go_num1)

					# If we get a match
					if (len(go_flags) > 0):
						go_rows+=1
						writer.writerow(row)

						for go_num in go_flags:
							for elem in elems:
								if elem == "taxid":
									taxids_conf = row["multi_taxids_confidence"].split(",")
									tc_list = []
									maximum = 0.0
									for taxid_conf in taxids_conf:
										tc_term = taxid_conf.split(":")
										if (float(tc_term[1]) > maximum):
											maximum = float(tc_term[1])
											tc_list.append(tc_term)
									for tc in tc_list:
										if float(tc[1]) >= maximum:
											go_dict[go_num][elem][tc[0]]+=1
								else:
									go_dict[go_num][elem][row[elem]]+=1


						for i in range(9):
							if (int(row[bpoc_names[i]]) > 0):
								counts[i]+=1
								for elem in elems:
									bpoc_dicts[i][elem][row[elem]]+=1

						


	# write to the text file
	txt_file = open(txtfilename, "w")
	txt_file.write("Total rows: " + str(total_rows) + "\n")
	txt_file.write("Rows with at least one of the GO numbers: " + str(go_rows) + "\n")
	txt_file.write("Percentage of go rows: " + str((1.0*go_rows)/total_rows) + "\n")

	txt_file.write("\n\n")



	# Go number information
	for go_num in go_dict:
		txt_file.write("Go num: " + go_num + "\n")
		for elem in go_dict[go_num]:
			txt_file.write(elem + ":\n")
			for i in go_dict[go_num][elem]:
				txt_file.write("key: " + str(i) + ", value: " + str(go_dict[go_num][elem][i]))
				txt_file.write("\n")
			txt_file.write("\n")
		txt_file.write("\n")

	txt_file.write("\n")


	# BPOC-go information

	for i in range(len(counts)):
		txt_file.write("Percentage of " + bpoc_names[i] + ": " + str((1.0*counts[i])/total_rows) + "\n")

	txt_file.write("\n\n")
	
	for i in range(len(bpoc_dicts)):
		for elem in elems:
			txt_file.write(bpoc_names[i] + " " + elem + " ")
			for j in bpoc_dicts[i][elem]:
				txt_file.write("key: " + str(j) + ", value: " + str(bpoc_dicts[i][elem][j]))
				txt_file.write("\n")
			txt_file.write("\n")
		txt_file.write("\n")



# filepath = "/Users/winnieli/Documents/summer2020microbes/"


# # Must input numbers with the "GO:" in front of it
# go_numbers = ["GO:0046718", "GO:0016032"]

# def go_query_filepath(filepath, filename, destfilename, txtfilename):
# 	go_query(go_numbers, filepath+filename, filepath+destfilename, filepath+txtfilename)

# # go_query_filepath(filepath, "S01_trim25_fast_seqscreen_report.tsv", "S01_go_revised.tsv", "S01_go.txt")
# # go_query_filepath(filepath, 'S02_trim25_fast_seqscreen_report.tsv', 'S02_go_revised.tsv', 'S02_go.txt')
# # go_query_filepath(filepath, 'S03_trim25_fast_seqscreen_report.tsv', 'S03_go_revised.tsv', 'S03_go.txt')
# # go_query_filepath(filepath, 'S04_trim25_fast_seqscreen_report.tsv', 'S04_go_revised.tsv', 'S04_go.txt')
# # go_query_filepath(filepath, 'SRR10903401_combined_trim25_fast_seqscreen_report.tsv', 'SRR401_go_revised.tsv', 'SRR401_go.txt')
# # go_query_filepath(filepath, 'SRR10903402_combined_trim25_fast_seqscreen_report.tsv', 'SRR402_go_revised.tsv', 'SRR402_go.txt')
# # go_query_filepath(filepath, 'SRR10971381_combined_trim25_seqscreen_report.tsv', 'SRR381_go_revised.tsv', 'SRR381_go.txt')


