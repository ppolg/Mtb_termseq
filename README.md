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

#### Mtb_DEseq.R
DEseq analysis with DEseq2 on RNAseq and term-seq data.

#### Mtb_NUGA.R
Analysis of overlapping ORFs/uORFs, especially of 4-nt potential TeRe regulatory overlaps. Downstream analysis.
Figures 5B-F, Supplementary Figure X

#### Mtb_PS_distance.R
Analyse and plot distances between 5' monoP and term-seq peaks. What window to use to differentiate TS and PPS? How does this compare to random chance?
Figure 1C

#### Mtb_peaksfrequency.R
Additinal analysis on the "profiles" of TTS, as discussed in text (related to the data in Figure 2)

#### Mtb_perdistance.R
Distribution of TTS Scores and RT scores throughout TTSs. Additional data for subgroups of these terminators.
Figure 4A and 4D, Supplementary Figures XYZ

#### Mtb_plotdepth.R
Nucleotide occupancy heatmaps and aggregate plot (for legibility) around TTSs (called from term-seq data), based on RNAseq coverage. Shows changes in read-through between different timepoints after induction of Rho depletion. Also for various subsets of data.
Figures 4C and 4F, Supplementary Figures XYZ

#### Mtb_pvalues.R
Statistical significance analysis. Heatmaps visualising different q-value and TTS-Score and/or RT Score cutoffs.
Supplementary Figures XY



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

Alexandre D'Halluin, for the initial work performed on this project,

Terry Kipkorir and Fabian Blombach for their great insights,

Kristine Arnvig and Irilenia Nobeli for their supervision, input and guiding hand,

Zaynah Patel for _something positive_, surely? Anything at all?

Caffeine, for being a legal stimulant

The BBSRC and the LIDo programme for providing me funding for my PhD, and giving me the opportunity to learn so much!
