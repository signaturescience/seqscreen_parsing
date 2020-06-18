#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 13:28:04 2020

@author: Student
"""

import copy 
import pandas as pd

def parse_conf(cell, conf): 
    #splits contents of multiple_taxid cell and removes taxid's below the input
    #confidence level
    
    #converts string into dictionary
    taxids = eval( "{" + cell + "}")
    final_tids = copy.deepcopy(taxids)
        
    #removes taxids below conf 
    for key, val in taxids.items(): 
        if val < conf: 
            final_tids.pop(key)
                
    if final_tids == {}: 
        return None    
    else: 
        temp = repr(final_tids).replace("{", "")
        strdict = temp.replace("}", "")
        return strdict
    
def remove_tied(cell, thresh): 
    #removes all taxids of equally high confidence from a cell
    
    #converts string into dictionary
    taxids = eval( "{" + cell + "}")    
    rev = {}
    
    for key, value in taxids.items():
        rev.setdefault(value, set()).add(key)
    
    for key in rev.keys(): 
        if len(rev[key]) > thresh: 
            for tid in rev[key]: 
                taxids.pop(tid) 
    if taxids == {}: 
        return None 
    else: 
        temp = repr(taxids).replace("{", "")
        strdict = temp.replace("}", "")
        return strdict
    
def remove_equal(cell):   
    return remove_tied(cell, 1)

def parse_funcs(inputname, destname, func, attrs): 
    """
    inputname : str representing input file path. 
    destname : str representing output file path. 
    func : function to be used to parse data (inputs and outputs str)
    attrs : tuple containing other inputs needed for func.

    parses data based on one of the given functions, writes edited data
    """
    df = pd.read_csv(inputname, sep="\t")
    l = "multi_taxids_confidence"
    df[l] = df[l].apply(func, 1, args=attrs)
    df = df.dropna()
    df.to_csv(destname, sep="\t", index=False)
    
# def parse_taxid(inputname, destname, conf): 
#     #creates new file only containing taxids above the input confidence 
#     df = pd.read_csv(inputname, sep="\t")
#     l = "multi_taxids_confidence"
#     df[l] = df[l].apply(parse_conf, 1, args=(conf,))
#     df = df.dropna()
#     df.to_csv(destname, sep="\t", index=False)

# print(parse_taxid("/Users/Student/Desktop/testinput.txt", "/Users/Student/Desktop/hebele.txt", .8))

def count_taxids(inputname, destname): 
    df = pd.read_csv(inputname, sep="\t")
    l = "taxid"
    
    #count frequency of each taxid. first element does regular count, second does relative count
    inp = [df[l].value_counts(), df[l].value_counts(normalize=True)]
    
    #concatonate count dataframes into one, add labels & total row. 
    summary = pd.concat(inp, 1)
    summary.loc['total']= summary.sum(numeric_only=True, axis=0)
    summary = summary.reset_index()
    summary.columns = ["taxid", "count", "percentage"]
    summary.to_csv(destname, sep="\t", index=False)
    
def assume_human(inputname, destname): 
    df = pd.read_csv(inputname, sep="\t")
    l = "multi_taxids_confidence"
    df[df[l].str.contains("9606:")]
    df.loc[df[l].str.contains("9606:"), "taxid"] = 9606
    
    return df 

#print(assume_human("/Users/Student/Desktop/ooooo.txt", "bo"))
    

