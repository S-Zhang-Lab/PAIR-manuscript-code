#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_CellBC   # Specify job name

module load python
source /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/venv/bin/activate

DATA_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10

python 02_Cellbc_Whitelist_Mapping.py \
    --whitelist ${DATA_DIR}/code/PAIR_extraction/mRNA_barcodes_with_prefix.csv \
    --extracted_data ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_mismatch01_AdjR2.txt \
    --output ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.mapped.txt \
    --unmapped ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.unmapped.txt \
    --tolerance 2 \
    --k 10 \
    --log_file ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.mapping.log.txt

