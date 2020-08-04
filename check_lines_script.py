import subprocess
import os
import gzip
import shutil
from Bio import SeqIO
import pandas as pd

OSF_DOWNLOAD_DIR = "7nrd3/osfstorage/Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_Filtered_Dataset/"
LOCAL_DIR = ""
# You will have to change this to /home/dbs/SeqScreenDB
DATABASE_LOCATION = "/home/dbs/SeqScreenDB"
# You will have to change this as well
USERNAME = "wwl3@rice.edu"

for filename in os.listdir(OSF_DOWNLOAD_DIR):
    # print("DF", df)
    # print(filename)
    if filename.endswith(".tsv"): 
        file = os.path.join(OSF_DOWNLOAD_DIR, filename)
        filename_short = filename[:-4]
        if os.path.exists(file):
            try:
                num_lines = subprocess.check_output("grep -c \"@\" {}".format(file), shell=True)
                # print(filename_short, num_lines)
            except subprocess.CalledProcessError as e:
                print(filename_short, e.output)
                
            

