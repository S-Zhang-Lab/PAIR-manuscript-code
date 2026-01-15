#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_repeating       # Specify job name

module load python
source /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/venv/bin/activate

DATA_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/PAIR_TT_C12
OUT_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/PAIR_output

# full run
python check_PAIR_repeating_seq.py \
  --fastq ${DATA_DIR}/PAIR_TT_C12-SI_TT_C12_2333FLLT3_S3_L007_R2_001.fastq.gz \
  --n_reads 0 \
  --max_mismatch 1 \
  --out_prefix ${OUT_DIR}/PAIR_TT_C12_L007_R2_whitelist_mm1_all


