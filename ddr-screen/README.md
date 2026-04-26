# PAIR-DDR

Investigation of Transcriptome Interaction in DNA Repair Pathway by Programmable CRISPR Paired Sequencing

### 01 Background

#### *Experiment name*

Systematically assay gene-gene interaction in DNA damage response by PAIR-seq

#### *Organism*

Homo sapiens

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

293T PAIR cells (a stable cell line capable of expressing Cas13d and CRISPRa) are cultured in 10 cm dishes. The DDR library (comprising a total of 10,525 different PAIR RNAs, including NC) is packaged into lentivirus and used to infect 293T PAIR cells at a low MOI of 0.2. After 48 hours of infection, 10 µg/mL puromycin is added for drug selection, and the medium is replaced once the blank control cells have died. After 7 days of drug selection, the DDR-infected cells are split into three dishes, followed by transfection with SaCas9/sgRNA and a BFP DNA donor. 4 days post-transfection, flow cytometry is used to sort BFP+ cells, and different samples are harvested for library construction.

#### *Groups*

1.  Biological Repeat 1: DDR1 (unsorted), DDR2 (BFP sorted)

2.  Biological Repeat 2: DDR3 (unsorted), DDR4 (BFP sorted)

3.  DDR library infected 293T cells after Puro selection: DDR5

#### *Amplicon sequencing*

The **sorted** cells and **unsorted** cells and puro selected cells are harvested and genomic DNA is extracted individually. The three groups of genomic DNA are used for library construction with different i5/i7 index. The amplicon library structure can be found here: <https://benchling.com/s/seq-jUHvMUn4ehvBeo0N6s9f?m=slm-zhwS7H7MuqeoNZrmSJCt> and under ./Misc folder [PDF](./Misc/DDR_nestPCR_libraray-sequence.pdf) and [GB](./Misc/DDR_nestPCR_libraray.gb) file.

***Experiment date*** 11/15/2024

***Platform*** NovaSeq X Plus (PE150) 375G of theoretical data per lane

***RAW data and info sheet:*** see [`../data/README.md`](../data/README.md)
for raw FASTQ accessions and processed-data download links.

### 02. Fastq raw data processing and count matrix generation

Data quick check and a series of raw fastq processing steps can be found in [UTSW25_DDR.sh](Codes/Fastq_processing_code/UTSW25_DDR.sh). The raw sequencing reads are QC'ed using fastp. The python scripts used to extract different elements and mapping steps can be found under /Codes/Fastq_processing_code/. The final raw count_matrix table ([kmer20_HD1_count_matrix.csv](./Misc/kmer20_HD1_count_matrix.csv)) and [kmer20_HD1_count_matrix_for_mageck.txt](./Misc/kmer20_HD1_count_matrix_for_mageck.txt) are located under ./Misc

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
python 02_assemble_count_matrix.py --whitelist_file ${PROJECT_ROOT:-.}/Sub_info/Final_oligos_withID.csv \
                                   --output_file kmer20_HD1_count_matrix.csv \
                                   --input_files DDR1_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR2_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR3_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR4_Final_kmer_mapped_count_HD1_T50_kmer20.csv DDR5_Final_kmer_mapped_count_HD1_T50_kmer20.csv

# step 5 MAGeCK (primary statistical method)
mageck test -k kmer20_HD1_count_matrix_for_mageck.txt \
    -t DDR2,DDR4 \
    -c DDR1,DDR3 \
    --norm-method total \
    -n DDR_default
```

### 03. DDR analysis in R

#### Statistical analysis overview

The PAIR count matrix from this combinatorial CRISPR screen is highly sparse (many zero counts), which is expected given the large number of unique PAIRs (10,525) and the stochastic nature of lentiviral infection and sorting. This sparsity has direct implications for statistical testing, and the choice of primary method is critical.

**DESeq2** is the **primary statistical method** for differential abundance testing. Among available parametric tools, DESeq2's shrinkage-based dispersion estimation is the most robust for n=2 designs, partially compensating for low replication through empirical Bayes moderation of per-gene variance. Significant PAIRs are identified at padj < 0.01 and |log2FC| > 5 with a pre-filter removing PAIRs with ≥3 zero counts.

**edgeR** serves as a **secondary validation** using a GLM framework with batch correction (FDR < 0.05, |logFC| > 5). The FDR threshold is slightly relaxed relative to DESeq2 because the GLM batch correction with n=2 leaves minimal residual degrees of freedom, making edgeR's FDR estimates overly conservative at stricter cutoffs. Results are additionally filtered for concordant fold-change direction across both biological replicates.

**MAGeCK** (Li et al., *Genome Biology*, 2014) is included as a **reference** for completeness but is not the primary method. While MAGeCK was designed for CRISPR screen data, its FDR calibration is unreliable for this focused library: it calls ~83% of tested PAIRs as significant at FDR < 0.05, indicating that the negative binomial model and robust rank aggregation are overfitting to this dataset. Two factors contribute: (1) the library is a focused combinatorial panel (10,525 PAIRs) rather than a genome-wide library with extensive sgRNA diversity per gene, and (2) MAGeCK's gene-level aggregation is not meaningful here because each "gene" is actually a gene pair with limited sgRNA representation.

**Gene-level frequency analysis** (script 05) identifies individual genes that appear at unexpectedly high frequency among DESeq2-significant PAIRs. A hypergeometric test compares the observed frequency of each gene against the expected frequency given its representation in the full library, providing a formal null model rather than reporting raw counts alone.

#### Pipeline execution

The R analysis pipeline (scripts 02–05) can be run as a single pipeline:

```r
# From RStudio with PAIR_DDR.Rproj open:
source("Codes/00_run_all.R")

# Or from the command line:
Rscript Codes/00_run_all.R
```

**NOTE:** MAGeCK output is optional. If `DDR_default.sgrna_summary.txt` and `DDR_default.gene_summary.txt` are placed in `./Results/MAGeCK/`, the integration script will include MAGeCK in the cross-validation. If absent, the pipeline runs with DESeq2 and edgeR only.

#### R scripts

1.  [02_count_matrix_analysis.R](./Codes/02_count_matrix_analysis.R): Exploratory analysis of the PAIR count matrix. Performs CPM-like normalization, computes log fold changes between sorted and unsorted cells within each replicate, identifies PAIRs with extreme concordant fold changes, and reformats the count matrix for MAGeCK input. Thresholds here are arbitrary and used for visualization only.

2.  [03_DE-seq2.R](./Codes/03_DE-seq2.R): **[Primary]** DESeq2 differential abundance analysis with batch correction (design: `~ batch + condition`). PAIRs with ≥3 zero counts are pre-filtered. Significant PAIRs identified at padj < 0.01 and |log2FC| > 5. Exports a full results table (`DEseq2_all_results.csv`) for cross-validation in script 05. Generates a heatmap and volcano plot with gene pair labels.

3.  [04_EdgeR.R](./Codes/04_EdgeR.R): **[Secondary]** edgeR GLM analysis with batch correction. FDR < 0.05 (relaxed relative to DESeq2 due to loss of degrees of freedom from batch correction with n=2), |logFC| > 5. Results are additionally filtered for concordant direction of change across both biological replicates. Exports a full results table (`edgeR_all_results.csv`) for cross-validation in script 05. Generates a heatmap and volcano plot.

4.  [05_integration_analysis.R](./Codes/05_integration_analysis.R): **[Primary]** DESeq2-based hit integration and cross-validation. Loads DESeq2 significant PAIRs as the primary hit set, cross-validates against edgeR and MAGeCK (if available), and performs gene-level frequency analysis with a hypergeometric null model to identify genes with enrichment/depletion rates exceeding chance expectation given library composition.

5.  [06_NBN_focus_analysis.R](./Codes/06_NBN_focus_analysis.R): Targeted analysis of NBN (nibrin) gene interactions from DESeq2-significant PAIRs. Generates two main visualizations: (1) enrichment bias panel showing the fraction of NBN-containing PAIRs enriched in BFP+ cells with statistical assessment (binomial test, Fisher's exact test, permutation test correcting for post-hoc selection), and (2) partner waterfall plot ranking all NBN-containing PAIRs by log2FC with color coding by DDR pathway membership (NHEJ, HR, DSB sensing, other). Outputs include individual and combined PDF figures plus results tables for pathway analysis. Uses official HUGO gene nomenclature throughout.

6.  [07_Target_Network_Analysis.R](./Codes/07_Target_Network_Analysis.R): Comprehensive network analysis and pathway positioning of the 20 target genes in the dual modulation library. Integrates STRINGdb protein-protein interactions with custom DDR pathway classification and KEGG enrichment to show (1) pathway cartography (genes organized by canonical pathway roles), (2) functional STRING network with pathway color-coding and network topology metrics, and (3) KEGG pathway enrichment results. Calculates network centrality measures (degree, betweenness, eigenvector centrality) and clustering coefficients to identify hub genes and cross-pathway bridges. Generates publication-ready visualizations showing library coverage of DDR pathway space and functional interaction landscape.

#### Deprecated scripts

-   `05_common_gene_analysis.R`: Original DESeq2 ∩ edgeR intersection analysis. Replaced by `05_integration_analysis.R`.
-   `05_MAGeCK_integration.R`: MAGeCK-primary integration analysis. Replaced by `05_integration_analysis.R` after determining that MAGeCK's FDR calibration is unreliable for this focused library (calls ~83% of PAIRs as significant).
