#!/usr/bin/env python3
# -*- coding: utf-8 -*-

def parse_bpoc(filename, destname): 
    
    # splits file into lines
    file = open(filename, "r")
    lines = file.readlines()
    splitlines = []
    
    # splits lines into lists
    for line in lines: 
        splitlines.append(line.split("\t"))

    # initializes counts (list of lists)
    end = []   
    
    # initializes dictionaries to store bpoc counts
    line_count = 0    
    count = {}
    count[14] = 0 
    
    for ind in range(5, 16): 
        count[ind] = 0
    
    labels = ['-', 'adhesion', 'secretion', 'host_cell_death', 'antibiotic', 
              'invasion', 'evasion', 'cytotoxicity', 'degrade_ecm', 
              'disable_organ', 'bpoc', 'total']
    
    # searches for lines with BPoCs 
    for line in splitlines: 
        if line[5] != "-": 
            count[15] += 1 
            line_count += 1 
            
        ind = 5
        
        # assumes 1 bpoc / line 
        while ind <= 13: 
            if line[ind] == "1" and line not in end: 
                end.append(line)
                count[ind] += 1
                count[14] += 1 
                #ind = 14
                ind += 1
            
            elif line[ind] == "1" and line in end: 
                count[ind] += 1
                ind += 1
            
            else: 
                ind += 1
    
    end.insert(0, splitlines[0])
    count[15] -= 1 
    
    # calculates percentages, adds into per   
    per = {}
    for bpoc in count.keys(): 
        per[bpoc] = count[bpoc] * 100 / (line_count -1) 
    
    # write data into output file 
    output = open(destname, 'w') 
    
    sum_labels = "-" + "\t".join(labels)
    
    scount = sorted(list(count.items()), key = lambda x: x[0])
    wline = ["count"]
    for ind, c in scount: 
        wline.append(str(c))
        
    pcount = sorted(list(per.items()), key = lambda x: x[0])
    pline = ["percentages"]
    for ind, c in pcount: 
        pline.append(str(c))
    
    output.writelines([sum_labels + "\n", "\t".join(wline) + "\n", "\t".join(pline)+ "\n", "-- \n"])
    
    write = []
    for l in end: 
        write.append('\t'.join(map(str, l)))
    output.writelines(write)
    
    #return count, per, end

#lists01 = get_bpoc("/Users/Student/opt/anaconda3/envs/microbes1/S02_trim25_fast_seqscreen_report.tsv")

# print(lists01)


#test = parse_bpoc(
#    "/Users/Student/opt/anaconda3/envs/microbes1/S01_trim25_fast_seqscreen_report.tsv",
#    "test2.txt")

# test2 = parse_bpoc("/Users/Student/Desktop/testinput.txt", "test3.txt")

# print(deneme2)