## reu2020_microbes

This repository contains a series of scripts for extracting data from the seqscreen output.
All files generated by this tool can be found in `outputs/` directory.

## BPoC Report

Taking in a .tsv file of SeqScreen output, creates:
- a revised.tsv file with only rows that have at least BPoC identified
- a summary.txt file with counts of BPoCs, percentages, and breakdowns of each BPoC
- a Krona Plot for the rows with BPoCs


```bash
python bpoc_parse.py <path_to_seqscreen_output>
```

**example:**
```bash
python bpoc_parse.py examples/SRR10903401_seqscreen_report.tsv
```

<br/><br/>


## GO Term Queries

Taking in a .tsv file of SeqScreen output and a list of GO terms, creates
these files for each GO term:
- a revised.tsv file with only rows that have that GO term
- a summary.txt file for that GO term
- a Krona Plot for that GO term

```bash
python go_term_parse.py <path_to_seqscreen_output> <go terms separated by spaces>
```

**example:**
```bash
python taxid_parse.py examples/SRR10903401_seqscreen_report.tsv 0016310 0016032 0003824
python taxid_parse.py examples/SRR10903401_seqscreen_report.tsv GO:0016310 GO:0016032 
```

<br/><br/>


## Taxon Abundances

Takes in a .tsv file of SeqScreen output. Creates the following files depending on which function is called: 
```bash
python taxid_parse.py <path_to_seqscreen_output> --function <function_attributes(ifany)>
```
- "count_taxid" creates an output file listing each taxid, the number of times it appears in the SeqScreen file, and its relative percentage compared to other taxids. 
```bash
python taxid_parse.py testinput.tsv --count_taxid
```
- "assume_human" assumes the sequence is human if the human taxi (9606) is encountered in the multi_taxid column. Returns a version of the input file with an edited “taxid” column.
```bash
python taxid_parse.py testinput.tsv --assume_human
```
- "all_tied" removes all taxids that share the same confidence level with at least one other taxid. Returns a version of the input file with the taxid and multi_taxid columns edited to show changes. Also returns a Krona plot based on the edited multi_taxid column.
```bash
python taxid_parse.py testinput.tsv --all_tied
```
- "thresh_tied" takes in an input threshold and removes taxids that share the same confidence level, if the number of tied taxids are above the input threshold. Returns a version of the input file with “taxids” and “multi_taxid” columns edited to show changes, and a Krona plot based on the edited multi_taxid column.
```bash
python taxid_parse.py testinput.tsv --thresh_tied 10
```
- "parse_conf" removes all taxids in the multi_taxid column whose confidences are below the input confidence level. Returns a version of the input file with “taxid” and “multi_taxid” columns edited to show removed taxids and a Krona plot based on the edited multi_taxid column.
```bash
python taxid_parse.py testinput.tsv --parse_conf 0.8
```
