# Mtb_termseq
Termseq additional analysis for  D'Halluin _et al._ (2023)

(Manuscript currently on BioRXiv -  https://doi.org/10.1101/2022.06.01.494293)

## Repository information

This repository contains the code used for a multiple steps of data analysis and figure generation (written in R and bash) for the manuscript "Term-seq reveals an abundance of conditional, Rho-dependent termination in _Mycobacterium tuberculosis_."  Code for the analysis was fully written by Peter Polgar.

### Input files
are in the folder NUGA_Data 

### Output folder
is NUGA_Out. These folder names are carry-over from the previous iteration of the work largely focusing on 4nt overlaps in _Mtb_. The folders contain **all** input/output.

### R folder
Contains the R scripts used to analyse and plot data. All scripts are annotated in-code.

#### Filename

### Bash folder
Contains bash scipts used to automate RNAseq and term-seq data processing:

#### calc_depth.sh
Uses samtools depth to calculate nucleotide-based coverage from bam files in a folder

#### calc_htseq.sh
Runs samtools index and htseq-count on all bam files in folder, for downstream DEseq analysis.

#### calc_nreads.sh
Runs samtools view to quickly chek number of reads for all bam files

#### split_strands.sh
Uses samtools view to split all bam files in folder into forward and reverse strands. 


## Acknowledgements

Alexandre D'Halluin, for the initial work performed on this project

Terry Kipkorir and Fabian Blombach for their great insights,

Kristine Arnvig and Irilenia Nobeli for their supervision, input and guiding hand,

Zaynah Patel for _something positive_, I'm sure?

Caffeine, for being a legal stimulant

The BBSRC and the LIDo programme for providing me funding for my PhD, and giving me the opportunity to learn so much!
