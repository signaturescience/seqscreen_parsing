import subprocess
import os

def generate_krona(input_file):
    #input_file = "SRR10903402_trim30_rnaspades_default_seqscreen_report.txt"

    input_prefix = input_file.split(".")[0]
    output_file = input_prefix + ".html"

    krona_input_list = []

    def sort_by_id(elem):
        return int(elem.split("/")[1].rstrip(".txt").split("_")[1].lstrip("DN"))

    with open(input_file, "r") as input:
        lines = input.readlines()
        if not os.path.exists(f"{input_prefix}_temp"):
            os.mkdir(f"{input_prefix}_temp")

        for line in lines[1:]:
            elements = line.split("\t")
            query_id = elements[0]
            multi_taxids = elements[3].split(",")
            with open(f"{input_prefix}_temp/{query_id}.txt", "w") as krona_input:
                for i in range(len(multi_taxids)):
                    q = str(i+1)
                    t = multi_taxids[i].split(":")[0]
                    try:
                        s = multi_taxids[i].split(":")[1]
                    except IndexError:
                        print(query_id, multi_taxids)
                        q = "0"
                        s = "0"

                    krona_input.write(q+"\t"+t+"\t"+s+"\n")

            krona_input_list.append(f"{input_prefix}_temp/{query_id}.txt")

    #krona_input_list.sort(key=sort_by_id)
    #print(krona_input_list)

    subprocess.run(["ktImportTaxonomy", "-q", "1", "-t", "2", "-s", "3"] + krona_input_list + ["-o", output_file])

# input_file_list = ["SRR10903402_trim30_rnaspades_fast_seqscreen_report.txt", "SRR10903402_trim30_megahit_fast_seqscreen_report.txt", "SRR10903402_trim30_trinity_default_seqscreen_report.txt", "SRR10903402_trim30_trinity_fast_seqscreen_report.txt", \
#     "SRR10903401_trim30_trinity_fast_seqscreen_report.txt", "SRR10903401_trim30_trinity_default_seqscreen_report.txt", "SRR10903401_trim30_rnaspades_fast_seqscreen_report.txt", "SRR10903401_trim30_megahit_fast_seqscreen_report.txt"]

input_file_list = ["S01_revised.tsv", "S02_revised.tsv", "S03_revised.tsv", "S04_revised.tsv", "SRR401_revised.tsv", "SRR402_revised.tsv", "SRR381_revised.tsv"]

for input_file in input_file_list:
    generate_krona(input_file)