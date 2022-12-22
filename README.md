# AlphaFold Data Base disorder trimming

This script removes long disordered segments from AlphaFold models. It uses [`AlphaFold-disorder`](https://github.com/BioComputingUP/AlphaFold-disorder) predictions as the criterion for trimming. It will return a trimmed copy of every PDB file in a collection, except for those predicted to be 100% disordered. 

Long disordered segments are defined as **>25 consecutive residues with a disorder score > 0.581** in either N->C or C->N sense. These default values derive from the default window size -- i.e. number of residues over which prediction is smoothed -- and the disorder threshold used by `AlphaFold-disorder`. 

## Usage

0. Run `AlphaFold-disorder` on your collection of models. 
1. Run 

`python trim_AFDB.py -i data/AlphaFoldDB/sample -p data/AlphaFold-disorder/sample_pred.tsv` 

## Comments

### External dependencies

`AlphaFold-disorder` requires DSSP (Define Secondary Structure of Proteins). You can install it from [here](https://ssbio.readthedocs.io/)

### Output PDB format

Trimmed PDBs are missing the rich headers from the original ADFB models. Instead, just the structure (with the original numbering) is returned, as written by [biopython's PDBIO](https://biopython.org/docs/1.75/api/Bio.PDB.PDBIO.html). 