# NLRclassify
Classification of NLRs based on an existing phylogenetic tree

These set of scripts use the Fast Approximate Tree Classification approach to classify proteins based on a preexisting phylogenetic tree. Initial scripts were ported from [https://github.com/shailen/FAT-CAT](https://github.com/shailen/FAT-CAT). Reference of the approach is the following:
> Afrasiabi, C., Samad, B., Dineen, D., Meacham, C. Sjölander, K., "The PhyloFacts FAT-CAT Webserver: Ortholog Identification and Function Prediction using Fast Approximate Tree Classification," ***Nucleic Acids Research*** 2013; [doi:10.1093/nar/gkt399](https://doi.org/10.1093/nar/gkt399)

## Scripts
*fat_cat_build.py* builds the initial hidden Markov model (HMM) library and runs a single input sequence against the library.

*fat_cat_map.py* runs a single sequence against a preexisting HMM library.

## Software and modules
The following software and modules are requires to run the suite of scripts:

1. Python 2.7
  * BioPython
2. HMMER 3.1b1

Other versions are likely to be functional, but the versions described above were used in the development of the scripts.

## Example
### Build initial HMM library
```
python fat_cat_build.py -i sample_alignment_IDs.txt -t sample_tree.ml -m sample_alignment.msa -s new_sequence.fasta  -c 8
```

### Run a sequence against a preexisting HMM library
```
python fat_cat_map.py -i sample_alignment_IDs.txt -s new_sequence.fasta -c 8
```
