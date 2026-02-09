#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe
#$ -pe smp 48 # reserving 48 cores
#$ -l h_vmem=300G
#$ -q long
#$ -N cellranger_mRNAc9_HTOc1

module load bio/cellranger/8.0.1

OUTDIR="/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/PAIRseq_Rd2/mRNAc9_HTOc1"
SAMPLE_ID="mRNAc9_HTOc1"
MULTI_CONFIG="/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/code/SeqAlignment_Rd2/mRNAc9_HTOc1/multi_config_c9_c1.csv"

# Go to the desired output directory
cd "$OUTDIR" || exit 1

cellranger multi \
  --id=${SAMPLE_ID}_Output \
  --csv=${MULTI_CONFIG} \
  --localcores=48 \
  --localmem=295

echo "✅ Cell Ranger multi completed at $(date)"