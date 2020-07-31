import subprocess
import os
import gzip
import shutil
from Bio import SeqIO

OSF_DOWNLOAD_DIR = "Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_Filtered_Dataset/"
OSF_UPLOAD_DIR = "Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_Biocuration_Summaries/"
LOCAL_DIR = ""
# You will have to change this to /home/dbs/SeqScreenDB
DATABASE_LOCATION = "/home/dbs/SeqScreenDB"
# You will have to change this as well
USERNAME = "wwl3@rice.edu"

# Change this to the files you want to run SeqScreen on
files = [
"SRR5787583_trim25_fast_combined_seqscreen_report.tsv",
"SRR5787584_trim25_seqscreen_report.tsv",
"SRR5787585_trim25_seqscreen_report.tsv",
"SRR5787586_trim25_fast_combined_seqscreen_report.tsv"
# "SRR11092058_filtered_1_reads_trim25_kraken2_unclassified_seqscreen_report.tsv"

]

for f in files:
    file = f[:-4]
    fname_remote = os.path.join(OSF_DOWNLOAD_DIR, f"{file}.tsv")
    fname = os.path.join(LOCAL_DIR, f"{file}.tsv")
    # subprocess.check_output("osf -u {} -p 7nrd3 fetch {} {}".format(
    #         USERNAME,
    #         fname_remote,
    #         fname),
    #     shell=True)
    
    # print(f"Running seqscreen on {file}:{count} samples")
    subprocess.check_output("python reu2020_microbes/bpoc_parse.py {}".format(
            fname),
        shell=True)
	#Upload back to OSF
    subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
            USERNAME,
            os.path.join("outputs", f"{file}_revised.tsv"),
            os.path.join(OSF_UPLOAD_DIR, f"{file}_revised.tsv")),
        shell=True)
    subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
            USERNAME,
            os.path.join("outputs", f"{file}_summary.txt"),
            os.path.join(OSF_UPLOAD_DIR, f"{file}_summary.txt")),
        shell=True)
    subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
            USERNAME,
            os.path.join("outputs", f"{file}_revised.tsv_krona.txt.html"),
            os.path.join(OSF_UPLOAD_DIR, f"{file}_revised.tsv_krona.txt.html")),
        shell=True)
