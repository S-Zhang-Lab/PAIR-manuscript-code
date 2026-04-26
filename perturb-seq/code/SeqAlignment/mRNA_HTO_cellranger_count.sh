#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe
#$ -pe smp 36 # reserving 36 cores
#$ -q long
#$ -N cellranger_mRNA_HTO

module load bio/cellranger/8.0.1

OUTDIR="/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10"
SAMPLE_ID="HTO_GEX_TT_A12"
MULTI_CONFIG="/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/code/SeqAlignment/multi_config.csv"

# Go to the desired output directory
cd "$OUTDIR" || exit 1

cellranger multi \
  --id=${SAMPLE_ID}_Output \
  --csv=${MULTI_CONFIG} \
  --localcores=36 \
  --localmem=220

echo "✅ Cell Ranger multi completed at $(date)"