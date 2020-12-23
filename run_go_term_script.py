import subprocess
import os
import gzip
import shutil
from Bio import SeqIO

OSF_DOWNLOAD_DIR = "Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_Filtered_Dataset/"
OSF_UPLOAD_DIR = "Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_GO_Term_Summaries/"
LOCAL_DIR = ""
# You will have to change this to /home/dbs/SeqScreenDB
DATABASE_LOCATION = "/home/dbs/SeqScreenDB"
# You will have to change this as well
USERNAME = "wwl3@rice.edu"

# Change this to the files you want to run SeqScreen on
files = [
"examples/SRR11059940_filtered_1_reads_trim25_seqscreen_report.tsv"
#    ,"SRR11059941_filtered_1_reads_trim25_seqscreen_report.tsv"
#    ,"SRR11059942_filtered_1_reads_trim25_seqscreen_report.tsv"
#    ,"SRR11059943_filtered_1_reads_trim25_seqscreen_report.tsv"
#    ,"SRR11059944_filtered_1_reads_trim25_1_seqscreen_report.tsv"
#    ,"SRR11059945_filtered_1_reads_trim25_1_seqscreen_report.tsv"
#    ,"SRR11059946_filtered_1_reads_trim25_1_seqscreen_report.tsv"
]

go_nums = []

for f in files:
    file = f[:-4]
    fname_remote = os.path.join(OSF_DOWNLOAD_DIR, f"{file}.tsv")
    fname = os.path.join(LOCAL_DIR, f"{file}.tsv")
    subprocess.check_output("osf -u {} -p 7nrd3 fetch {} {}".format(
            USERNAME,
            fname_remote,
            fname),
        shell=True)
    
    # print(f"Running seqscreen on {file}:{count} samples")
    subprocess.check_output("python reu2020_microbes/go_term_parse.py {} {}".format(
            fname,
            go_nums.join(" "),
        shell=True))
	#Upload back to OSF
    for go_num in go_nums:
        subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
                USERNAME,
                os.path.join("outputs", f"{file}_go_num_{go_num}_revised.tsv"),
                os.path.join(OSF_UPLOAD_DIR, f"{file}_go_num_{go_num}_revised.tsv")),
            shell=True)
        subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
                USERNAME,
                os.path.join("outputs", f"{file}_go_num_{go_num}_summary.txt"),
                os.path.join(OSF_UPLOAD_DIR, f"{file}_go_num_{go_num}_summary.txt")),
            shell=True)
        subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
                USERNAME,
                os.path.join("outputs", f"{file}_go_num_{go_num}_revised.tsv_krona.txt.html"),
                os.path.join(OSF_UPLOAD_DIR, f"{file}_go_num_{go_num}_revised.tsv_krona.txt.html")),
            shell=True)
