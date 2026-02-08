#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_f2       # Specify job name

module load python
source /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/venv/bin/activate

DATA_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/PAIRseq_Rd2/01.RawData/CC_PAIR_f2
OUT_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/PAIRseq_Rd2/PAIR_f2

# full run
python 01_Extract_PAIR.py   --input_r1 ${DATA_DIR}/CC_PAIR_f2-SI_TT_F2_23357TLT4_S6_L004_R1_001.fastq.gz \
                            --input_r2 ${DATA_DIR}/CC_PAIR_f2-SI_TT_F2_23357TLT4_S6_L004_R2_001.fastq.gz \
                            --output_r1 ${OUT_DIR}/PAIR_R1_f2_mismatch01_AdjR2.fastq.gz \
                            --output_r2 ${OUT_DIR}/PAIR_R2_f2_mismatch01_AdjR2.fastq.gz \
                            --output_table ${OUT_DIR}/PAIR_f2_table_mismatch01_AdjR2.txt \
                            --whitelist /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/code/PAIR_extraction_Rd2/PAIR_perturb_whitelist.csv

