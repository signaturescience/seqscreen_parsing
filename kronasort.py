#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 15:59:43 2020

@author: Student
"""

import pandas as pd 
import os 

def krona_input(infile, outfile): 
    temp = pd.read_csv(infile, sep="\t")
    df = temp["multi_taxids_confidence"].str.split(",")
    df = pd.concat([temp["query"], df], 1)
    df = df.explode("multi_taxids_confidence")
    final = pd.concat([df["query"], df["multi_taxids_confidence"].str.split(":", expand=True)], 1)
    final.to_csv(outfile, sep="\t", header=False, index=False)
    
