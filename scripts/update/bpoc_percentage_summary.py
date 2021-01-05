import subprocess
import os
import gzip
import shutil
from Bio import SeqIO
import pandas as pd

OSF_DOWNLOAD_DIR = "7nrd3/osfstorage/Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_Biocuration_Summaries/"
OSF_UPLOAD_DIR = "Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_Biocuration_Summaries/"
LOCAL_DIR = ""
# You will have to change this to /home/dbs/SeqScreenDB
DATABASE_LOCATION = "/home/dbs/SeqScreenDB"
# You will have to change this as well
USERNAME = "wwl3@rice.edu"

cols = ["sample", "adhesion", "secretion", "host_cell_death",
             "antibiotic", "invasion", "evasion",
             "cytotoxicity", "degrade_ecm", "disable_organ"]
df = pd.DataFrame(columns = cols)
# df.set_index(["sample", "adhesion", "secretion", "host_cell_death",
#              "antibiotic", "invasion", "evasion",
#              "cytotoxicity", "degrade_ecm", "disable_organ"])

# subprocess.check_output("osf -u {} -p 7nrd3 fetch {} {}".format(
#         USERNAME,
#         os.path.join(OSF_DOWNLOAD_DIR, "bpoc_percentage_summary.tsv"),
#         os.path.join(OSF_UPLOAD_DIR, "bpoc_percentage_summary.tsv")),
#     shell=True)


for filename in os.listdir(OSF_DOWNLOAD_DIR):
    # print("DF", df)
    # print(filename)
    if filename.endswith(".txt"): 
        file = os.path.join(OSF_DOWNLOAD_DIR, filename)
        filename_short = filename[:-4]
        values = {"sample": filename_short}
        if os.path.exists(file):
            f = open(file, "r")
            f.readline()
            for i in range(1,10):
                data = f.readline().split()
                percentage = data[2]
                values[cols[i]] = percentage
            # df_values = pd.DataFrame.from_dict(values)
            df = df.append(values, ignore_index=True)
            f.close()

df.to_csv(os.path.join("outputs", "bpoc_percentage_summary.tsv"),
         sep='\t', index=False)
 
subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
        USERNAME,
        os.path.join("outputs", "bpoc_percentage_summary.tsv"),
        os.path.join(OSF_UPLOAD_DIR, "bpoc_percentage_summary.tsv")),
    shell=True)
