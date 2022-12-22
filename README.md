# AlphaFold Data Base disorder trim

This script removes long disordered segments from AlphaFold models. It uses [AlphaFold-disorder](https://github.com/BioComputingUP/AlphaFold-disorder) predictions as the criterion for trimming.

Long disordered segments are defined as >25 consecutive residues with a disorder score > 0.581. These values derive from the default window size -- i.e. number of residues over which prediction is smoothed -- and the disorder threshold used to calculate relative solvent accessibility (RSA) in AlphaFold-disorder. The rolling window is considered in both N->C and C->N senses.

## Usage