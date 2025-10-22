#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe
#$ -pe smp 24 #reserving 24 cores
#$ -q long
#$ -N HTO_lib

module load conda

conda activate citeseq


CITE-seq-Count \
  -R1 /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/HTO_TN_A3/HTO_TN_A3-SI_TN_A3_2333FLLT3_S1_L007_R1_001.fastq.gz \
  -R2 /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/HTO_TN_A3/HTO_TN_A3-SI_TN_A3_2333FLLT3_S1_L007_R2_001.fastq.gz \
  -t  /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/code/SeqAlignment_Separate/HTO_tags.csv \
  -cbf 1 \
  -cbl 16 \
  -umif 17 \
  -umil 28 \
  -cells 40000 \
  -o /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/HTO_Output