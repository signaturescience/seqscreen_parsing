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

	names = ["adhesion", "secretion", "host_cell_death",
			"antibiotic", "invasion", "evasion",
			"cytotoxicity", "degrade_ecm", "disable_organ"]

	# initialize our dictionaries
	adhesion_dict = {}
	secretion_dict = {}
	host_cell_death_dict = {}
	antibiotic_dict = {}
	invasion_dict = {}
	evasion_dict = {}
	cytotoxicity_dict = {}
	degrade_ecm_dict = {}
	disable_organ_dict = {}

	dicts = [adhesion_dict, secretion_dict, host_cell_death_dict, 
			antibiotic_dict, invasion_dict, evasion_dict,
			cytotoxicity_dict, degrade_ecm_dict, disable_organ_dict]

	elems = ["taxid", "organism", "gene_name", "uniprot", "uniprot_evalue"]

	for this_dict in dicts:
		for elem in elems:
			this_dict[elem] = defaultdict(int)

	# read and parse the file
	file = open(filename, "r")
	info = file.readlines()

	dest_file = open(destfilename, "w")

	for line in info:
		data = (line.split('\t'))
		data2 = data[2].split(";")
		go_flag = False
		for elem in data2:
			for go_num in go_numbers:
				if (elem[3:] == go_num):
					go_flag = True
		if data[5]!='-' and data[0]!="query":
			total_rows+=1
		if (go_flag):
			for i in range(9):
				if (int(data[i+5]) > 0):
					counts[i]+=1
					dicts[i]["taxid"][data[1]]+=1
					dicts[i]["organism"][data[15]]+=1
					dicts[i]["gene_name"][data[16]]+=1
					dicts[i]["uniprot"][data[17]]+=1
					dicts[i]["uniprot_evalue"][data[18]]+=1
			go_rows+=1
			dest_file.write(line)


	# write to the text file
	txt_file = open(txtfilename, "w")
	txt_file.write("Total rows: " + str(total_rows) + "\n")
	txt_file.write("Rows with at least one of the GO numbers: " + str(go_rows) + "\n")
	txt_file.write("Percentage of go rows: " + str((1.0*go_rows)/total_rows) + "\n")

	txt_file.write("\n\n")

	for i in range(len(counts)):
		txt_file.write("Percentage of " + names[i] + ": " + str((1.0*counts[i])/total_rows) + "\n")

	txt_file.write("\n\n")
	
	for i in range(len(dicts)):
		for elem in elems:
			txt_file.write(names[i] + " " + elem + " ")
			for j in dicts[i][elem]:
				txt_file.write("key: " + str(j) + ", value: " + str(dicts[i][elem][j]))
				txt_file.write("\n")
			txt_file.write("\n")
		txt_file.write("\n")


filepath = "/Users/winnieli/Documents/summer2020microbes/"
go_numbers = ["0055114", "0016491"]
go_query(go_numbers, filepath + "S01_trim25_fast_seqscreen_report.tsv", filepath + "S01_go_revised.tsv", filepath + "S01_go.txt")
# go_query('S02_trim25_fast_seqscreen_report.tsv', 'S02_revised.tsv', 'S02.txt')
# go_query('S03_trim25_fast_seqscreen_report.tsv', 'S03_revised.tsv', 'S03.txt')
# go_query('S04_trim25_fast_seqscreen_report.tsv', 'S04_revised.tsv', 'S04.txt')
# go_query('SRR10903401_combined_trim25_fast_seqscreen_report.tsv', 'SRR401_revised.tsv', 'SRR401.txt')
# go_query('SRR10903402_combined_trim25_fast_seqscreen_report.tsv', 'SRR402_revised.tsv', 'SRR402.txt')
# go_query('SRR10971381_combined_trim25_seqscreen_report.tsv', 'SRR381_revised.tsv', 'SRR381.txt')


