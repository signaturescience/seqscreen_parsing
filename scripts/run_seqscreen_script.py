import subprocess
import os
import gzip
import shutil
from Bio import SeqIO

OSF_DOWNLOAD_DIR = "Human_Metatranscriptomes/Final_Results/Pre-Processing/Trimmed_Filtered_Reads/"
OSF_UPLOAD_DIR = "Human_Metatranscriptomes/Final_Results/SeqScreen/SeqScreen_Filtered_Dataset/"
LOCAL_DIR = ""
# You will have to change this to /home/dbs/SeqScreenDB
DATABASE_LOCATION = "/home/dbs/SeqScreenDB"
# You will have to change this as well
USERNAME = "wwl3@rice.edu"

files = ["CRR125952_filtered_1_reads_trim25_1.fq.gz",
"CRR125953_filtered_1_reads_trim25_1.fq.gz",
"CRR125954_filtered_1_reads_trim25_1.fq.gz",
"CRR125955_filtered_1_reads_trim25_1.fq.gz",
"CRR125956_filtered_1_reads_trim25_1.fq.gz",
"CRR125957_filtered_1_reads_trim25_1.fq.gz",
"CRR125958_filtered_1_reads_trim25_1.fq.gz",
"CRR125959_filtered_1_reads_trim25_1.fq.gz"]

# i represents the sample number. For the S*_filtered datasets, its easy since i=[2,22],
# but for other filename sets it might be a bit harder
for f in files:
    file = f[:-6]
    fname_fq_compressed_remote = os.path.join(OSF_DOWNLOAD_DIR, f"{file}.fq.gz")
    fname_fq_compressed = os.path.join(LOCAL_DIR, f"{file}.fq.gz")
    fname_fq = os.path.join(LOCAL_DIR, f"{file}.fq")
    fname_fa = os.path.join(LOCAL_DIR, f"{file}.fa")
    subprocess.check_output("osf -u {} -p 7nrd3 fetch {} {}".format(
            USERNAME,
            fname_fq_compressed_remote,
            fname_fq_compressed),
        shell=True)
    # gunzip
    with gzip.open(fname_fq_compressed, 'rb') as f_in:
        with open(fname_fq, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    # convert from fastq to fasta. Probably faster to do via sed...
    with open(fname_fq) as f_in:  
        with open(fname_fa, 'w') as f_out:
            sequences = SeqIO.parse(f_in, "fastq")
            count = SeqIO.write(sequences, f_out, "fasta")
    # print(f"Running seqscreen on S{i:02d}:{count} samples")
    subprocess.check_output("seqscreen_fast.nf --fasta {} --databases {} --working {} --threads 60".format(
            fname_fa,
            DATABASE_LOCATION,
            file),
        shell=True)
	#Upload back to OSF
    subprocess.check_output("osf -u {} -p 7nrd3 upload {} {}".format(
            USERNAME,
            os.path.join("output_dir", "report_generation", "seqscreen_report.tsv"),
            os.path.join(OSF_UPLOAD_DIR, f"{file}_seqscreen_report.tsv")),
        shell=True)
