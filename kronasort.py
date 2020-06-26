#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:59:43 2020

@author: Student
"""

import pandas as pd 
import subprocess

def make_krona(infile): 
    temp = pd.read_csv(infile, sep="\t")
    df = temp["multi_taxids_confidence"].str.split(",")
    df = pd.concat([temp["query"], df], 1)
    df = df.explode("multi_taxids_confidence")
    final = pd.concat([df["query"], df["multi_taxids_confidence"].str.split(":", expand=True)], 1)
    input_prefix = infile.split(".")[0]
    tsvname = f"{input_prefix}_kt_in.tsv"
    final.to_csv(tsvname, sep="\t", header=False, index=False)
    subprocess.run(["ktImportTaxonomy", "-q", "1", "-t", "2", "-s", "3", tsvname, "-o", f"{input_prefix}krona.html"])
    
