## reu2020_microbes

This repository contains a series of scripts for extracting data from the seqscreen output.
All files generated by this tool can be found in `outputs/` directory.

## Installation and Dependencies

Users must install [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html) 
and [krona](https://github.com/marbl/Krona/wiki/Installing).


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
python go_term_parse.py examples/SRR10903401_seqscreen_report.tsv 0016310 0016032 0003824
python go_term_parse.py examples/SRR10903401_seqscreen_report.tsv GO:0016310 GO:0016032 
```

<br/><br/>


## Taxon Abundances

Takes in a .tsv file of SeqScreen output. Creates the following files depending on which function is called: 
```bash
python taxid_parse.py <path_to_seqscreen_output> [--parse_conf PARSE_CONF][--thresh_tied THRESH_TIED] [--all_tied][--assume_human] [--count_taxid]
```
- "count_taxid" creates an output file listing each taxid, the number of times it appears in the SeqScreen file, and its relative percentage compared to other taxids. 
```bash
python taxid_parse.py examples/testinput.tsv --count_taxid
```
- "assume_human" assumes the sequence is human if the human taxi (9606) is encountered in the multi_taxid column. Returns a version of the input file with an edited “taxid” column.
```bash
python taxid_parse.py examples/testinput.tsv --assume_human
```
- "all_tied" removes all taxids that share the same confidence level with at least one other taxid. Returns a version of the input file with the taxid and multi_taxid columns edited to show changes. Deletes rows where all the taxids in the multi_taxid column are removed by the function.
```bash
python taxid_parse.py examples/testinput.tsv --all_tied
```
- "thresh_tied" takes in an input threshold and removes taxids that share the same confidence level, if the number of tied taxids are above the input threshold. Returns a version of the input file with “taxids” and “multi_taxid” columns edited to show changes, and a Krona plot based on the edited multi_taxid column. Deletes rows where all the taxids in the multi_taxid column are removed by the function.
```bash
python taxid_parse.py examples/testinput.tsv --thresh_tied 2
```
- "parse_conf" removes all taxids in the multi_taxid column whose confidences are below the input confidence level. Returns a version of the input file with “taxid” and “multi_taxid” columns edited to show removed taxids and a Krona plot based on the edited multi_taxid column. Deletes rows where all the taxids in the multi_taxid column are removed by the function.
```bash
python taxid_parse.py examples/testinput.tsv --parse_conf 0.8
```

## Report Visualisation 

Contains visualisation tools for use on SeqScreen outputs and Kraken2 Reports, showing the distribution of taxa across multiple samples. 

Users must install [altair](https://altair-viz.github.io/getting_started/installation.html) in addition to the dependencies listed above.

## Visualizing Kraken2 Reports

Takes in the following and creates an interactive stacked, grouped bar chart. 
1. the path to a directory containing only Kraken2 reports (except for hidden files)
2. the path to a file containingg grouping information (see sample for formatting)
3. taxids, separated by one space.


```bash
python report_vis.py <path_to_directory> -k  <path_to_groups> <taxids>
```

**example:**
```bash
python report_vis.py examples/Filtered_Kraken2_Reports -k examples/Sample_Groupings.txt 9606 10847
```

<br/><br/>


## Visualising SeqScreen Outputs 

Taking in the path to a directory containing only SeqScreen outputs (except for hidden files), creates an interactive stacked bar chart. By default, it visualises taxa using taxonomic ranks assigned by SeqScreen.


```bash
python report_vis.py <path_to_directory> -s
```

Takes optional arguments (specifying taxonomic rank and email address), allowing user to specify the taxonomic rank of organisms in the report. Ranks in SeqScreen outputs can vary, but usually identify strains. This is a useful flag to use for visualisation on the species level.

The specified string must be a valid string in the NCBI taxonomic database. Email address is required to use Entrez, which is used to fetch taxonomic information.

**example:**
```bash
python report_vis.py examples/seqscreen_outputs -s <rank> <email_address>
```

```bash
python report_vis.py examples/seqscreen_outputs -s species your@address.com
```

