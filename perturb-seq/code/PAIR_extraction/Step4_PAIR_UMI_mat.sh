#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_UMI   # Specify job name

module load python
source /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/venv/bin/activate

DATA_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10

python 04_Generate_PAIR_UMI_matrix.py \
    -i ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped.txt \
    -s ${DATA_DIR}/PAIR_output/PAIR_matched_sparse_UMI_matrix \
    -f ${DATA_DIR}/PAIR_output/PAIR_matched_frequency.txt
