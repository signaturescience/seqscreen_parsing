# reu2020_microbes
## Taxon Abundances

Takes in a .tsv file of SeqScreen output. Creates the following files depending on which function is called: 

- "count_taxid" creates an output file listing each taxid, the number of times it appears in the SeqScreen file, and its relative percentage compared to other taxids. 
- "assume_human" assumes the sequence is human if the human taxi (9606) is encountered in the multi_taxid column. Returns a version of the input file with an edited “taxid” column.
- "all_tied" removes all taxids that share the same confidence level with at least one other taxid. Returns a version of the input file with the taxid and multi_taxid columns edited to show changes. Also returns a Krona plot based on the edited multi_taxid column.
- "multi_tied" takes in an input threshold and removes taxids that share the same confidence level, if the number of tied taxids are above the input threshold. Returns a version of the input file with “taxids” and “multi_taxid” columns edited to show changes, and a Krona plot based on the edited multi_taxid column.
- "parse_conf" removes all taxids in the multi_taxid column whose confidences are below the input confidence level. Returns a version of the input file with “taxid” and “multi_taxid” columns edited to show removed taxids and a Krona plot based on the edited multi_taxid column.

**example:**
```bash
python taxidparse.py <path_to_seqscreeen_output> <function_name> <function_attributes (if any)> 
```
**additional flags for some functions:**

  - "multi_tied" takes in the maximum number of tied taxids that share the same confidence allowed.
  ```bash
    python taxidparse.py <path_to_seqscreeen_output> multi_tied -t <threshold>
```
  - "parse_conf" takes in the cutoff confidence level for taxids in the "multi_taxid" column.
  ```bash
    python taxidparse.py <path_to_seqscreeen_output> parse_conf -c <confidence>
```
