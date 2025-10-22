#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe
#$ -pe smp 24 #reserving 24 cores
#$ -q long
#$ -N mRNA_lib

module load bio/cellranger/8.0.1

OUTDIR="/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10"
SAMPLE_ID=mRNA_LIB
fastq_path=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/mRNA_TT_A12
Transcriptome=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/ref_genome/refdata-gex-GRCh38-2024-A

# Go to the desired output directory
cd "$OUTDIR" || exit 1

cellranger count --id=${SAMPLE_ID}_Output \
 --fastqs=${fastq_path} \
 --transcriptome=${Transcriptome} \
 --localcores=24 \
 --create-bam=false
