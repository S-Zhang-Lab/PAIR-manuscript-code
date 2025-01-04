# PAIR-DDR

Investigation of Transcriptome Interaction in DNA Repair Pathway by Programmable CRISPR Paired Sequencing

### Background introduction (CC add)

U6-driven PAIR RNA cassette structure can be find here: <https://benchling.com/s/seq-SQ0Ur1pqXan8cwLW510g?m=slm-h2zfzbW6z9kpwVhFuXMc>

### Design of Dual crRNA cassette

The dual regulatory crRNA cassette is with the following structure:

e1 \<- "AGGGCCTATTTCCCATGATTcgtctcacaccg"

e2 \<- "N20": CRISPRa gRNA

e3 \<- "gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAAC"

e4 \<- "n23": Cas13d (CasRx) crRNA

e5 \<- "TTTTTTT"

e6 \<- "BC14": 14 digit BC with hamming distance more than 5, which allow 2 bit error correction

e7 \<- "ctacagagacgcacttgtacttcagcggtc"

### Experimental details

293T PAIR cells (a stable cell line capable of expressing Cas13d and CRISPRa) are cultured in 10 cm dishes. The DDR library (comprising a total of 10,525 different PAIR RNAs, including NC) is packaged into lentivirus and used to infect 293T PAIR cells at a low MOI of 0.2. After 48 hours of infection, 10 µg/mL puromycin is added for drug selection, and the medium is replaced once the blank control cells have died. After 7 days of drug selection, the DDR-infected cells are split into three dishes, followed by transfection with SaCas9/sgRNA and a BFP DNA donor. 4 days post-transfection, flow cytometry is used to sort BFP+ cells, and different samples are harvested for library construction.

### Groups

Biological Repeat 1: DDR1 (unsorted), DDR2 (BFP sorted) Biological Repeat 2: DDR3 (unsorted), DDR4 (BFP sorted) DDR library infected 293T cells after Puro selection: DDR5

### Amplicon sequencing

The sorted cells and unsorted cells and puro selected cells are harvested. Then extract genomic DNA individually. The three groups of genomics DNA are using for library construction with different i5/i7 index. The amplicaon library structure can be found here: <https://benchling.com/s/seq-jUHvMUn4ehvBeo0N6s9f?m=slm-zhwS7H7MuqeoNZrmSJCt>

### Experiment date and platform

11/15/2024, Platform NovaSeq X Plus (PE150)

### Fastq data processing and count matrix generation

The raw sequencing reads QC'ed using fastp. The fastq processing steps can be found in [UTSW25_DDR.sh](Codes/Fastq_processing_code/UTSW25_DDR.sh) The python scripts used to extract different elements and mapping steps can be found under /Codes/Fastq_processing_code/ The final raw count_matrix table is located at /Misc/count_matrix.csv

### Preliminary analysis in R

Count normalization, manual DE gene calculation and logFC cutoff etc. has be documented in this [R script](Codes/02_count_matrix_analysis.R). The exported count matrix is [here](Results/PAIR_DDR_count_matrix.csv).
