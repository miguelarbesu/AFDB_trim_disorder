# AlphaFold Data Base disorder trim

This script trims long disordered segments from AlphaFold models. It uses [`AlphaFold-disorder`](https://github.com/BioComputingUP/AlphaFold-disorder) predictions as the criterion for trimming.

Long disordered segments are defined as **>25 consecutive residues with a disorder score > 0.581** in either N->C or C->N sense. These default values derive from the default window size -- i.e. number of residues over which prediction is smoothed -- and the disorder threshold used by `AlphaFold-disorder`. 

## Usage

`python trim_AFDB.py -i data/AlphaFoldDB/sample -p data/AlphaFold-disorder/sample_pred.tsv` 
