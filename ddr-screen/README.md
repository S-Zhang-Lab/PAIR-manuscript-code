# PAIR-DDR

Investigation of Transcriptome Interaction in DNA Repair Pathway by Programmable CRISPR Paired Sequencing

### 01 Backgound

#### *Experiment name*

Systematically assay gene-gene interaction in DNA damage response by PAIR-seq

#### *Organism*

Homo sapiens Experiment type Other

#### *Design of Dual crRNA cassette*

U6-driven PAIR RNA cassette structure can be find here: <https://benchling.com/s/seq-SQ0Ur1pqXan8cwLW510g?m=slm-h2zfzbW6z9kpwVhFuXMc>

The dual regulatory crRNA cassette is with the following structure:

1.  e1 \<- "AGGGCCTATTTCCCATGATTcgtctcacaccg"

2.  e2 \<- "N20": CRISPRa gRNA

3.  e3 \<- "gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAAC"

4.  e4 \<- "n23": Cas13d (CasRx) crRNA

5.  e5 \<- "TTTTTTT"

6.  e6 \<- "BC14": 14 digit BC with hamming distance more than 5, which allow 2 bit error correction

7.  e7 \<- "ctacagagacgcacttgtacttcagcggtc"

The DNA oligo library was designed using [01_Lib_assemble.R](./Codes/01_Lib_assemble.R)

#### *Experiment summary*

293T PAIR cells (a stable cell line capable of expressing Cas13d and CRISPRa) are cultured in 10 cm dishes. The DDR library (comprising a total of 10,525 different PAIR RNAs, including NC) is packaged into lentivirus and used to infect 293T PAIR cells at a low MOI of 0.2. After 48 hours of infection, 10 µg/mL puromycin is added for drug selection, and the medium is replaced once the blank control cells have died. After 7 days of drug selection, the DDR-infected cells are split into three dishes, followed by transfection with SaCas9/sgRNA and a BFP DNA donor. 4 days post-transfection, flow cytometry is used to sort BFP+ cells, and different samples are harvested for library construction. "

#### *Groups*

Biological Repeat 1: DDR1 (unsorted), DDR2 (BFP sorted) Biological Repeat 2: DDR3 (unsorted), DDR4 (BFP sorted) DDR library infected 293T cells after Puro selection: DDR5

#### *Amplicon sequencing*

The **sorted** cells and **unsorted** cells and puro seleted cells are harvest and extract genomic DNA individually. The three groups of genomics DNA are using for library construction with different i5/i7 index.The amplicon library structure can be found here: <https://benchling.com/s/seq-jUHvMUn4ehvBeo0N6s9f?m=slm-zhwS7H7MuqeoNZrmSJCt> and under ./Misc folder [PDF](./Misc/DDR_nestPCR_libraray-sequence.pdf) and [GB](./Misc/DDR_nestPCR_libraray.gb) file.

***Experiment date*** 11/15/2024

***Platform*** NovaSeq X Plus (PE150) 375G of theoretical data per lane

***RAW data and info sheet:*** RAW data is located on S.Zhang.Box shared drive:

1.  Info sheets: <https://utsw.box.com/s/redff0ap1bofn6ifszhnrynp5pmxme56>

2.  RAW data: <https://utsw.box.com/s/fd1l03t5hai81x4haeph3hlfwmgzgzhf>

### 02. Fastq raw data processing and count matrix generation

Data quick check and a series of raw fastq processing steps can be found in [UTSW25_DDR.sh](Codes/Fastq_processing_code/UTSW25_DDR.sh) The raw sequencing reads QC'ed using fastp. The python scripts used to extract different elements and mapping steps can be found under /Codes/Fastq_processing_code/ The final raw count_matrix table ([kmer20_HD1_count_matrix.csv](./Misc/kmer20_HD1_count_matrix.csv)) and [kmer20_HD1_count_matrix_for_mageck.txt](./Misc/kmer20_HD1_count_matrix_for_mageck.txt) are located under ./Misc

```{bash}
# step 1: extract DDR and BC
python UTSW25_DDR_BC_extraction.py -r1 DDR1_CKDL240042575-1A_22TFMKLT3_L7_1.fq.gz -r2 DDR1_CKDL240042575-1A_22TFMKLT3_L7_2.fq.gz -o DDR1_extraction_output.txt -l DDR1_extraction_output.log -p 100000

# step 2 calculate_frequencies
python calculate_frequencies.py -i DDR1_extraction_output.txt -o DDR1_unique_combinations_sorted.txt -v

# step 3 mapping (with k-mer error correction) to DDR gRNA whitelist and generate matrix
python 01_Mapping.py --input_file DDR1_unique_combinations_sorted.txt \
                     --whitelist_file Whitelist_withID.csv \
                     --output_file DDR1_Final_kmer_mapped_count_HD1_T50_kmer20.csv \
                     --tolerance 1 \
                     --count_threshold 50
                     
# step 4 assemble final count matrix
python 02_assemble_count_matrix.py --whitelist_file /project/pathology/SiZhang_lab/shared/Active_Projects/UTSW25_CC_DDR/Sub_info/Final_oligos_withID.csv \
                                   --output_file kmer20_HD1_count_matrix.csv \
                                   --input_files DDR1_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR2_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR3_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR4_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR5_Final_kmer_mapped_count_HD1_T50_kmer20.csv

# step 5 quick test using using MAGack. NOTE, the output file kmer20_HD1_count_matrix.csv should be formatted as the following based on Magack count table format. 
mageck test -k kmer20_HD1_count_matrix_for_mageck.txt \
    -t DDR2,DDR4 \
    -c DDR1,DDR3 \
    --norm-method total \
    -n DDR_default 
```

### 03. DDR analysis in R

{need to update}
