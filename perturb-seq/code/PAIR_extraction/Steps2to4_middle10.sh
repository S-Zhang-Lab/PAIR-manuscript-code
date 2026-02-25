#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_Steps2to4_middle10   # Specify job name

set -e  # stop on errors

module load python
source /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/venv/bin/activate

DATA_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10

python 02_Cellbc_Whitelist_Mapping.py \
    --whitelist ${DATA_DIR}/code/PAIR_extraction/mRNA_barcodes_with_prefix.csv \
    --extracted_data ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_mismatch01_AdjR2_middle10.txt \
    --output ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_CellBC.mapped_middle10.txt \
    --unmapped ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_CellBC.unmapped_middle10.txt \
    --tolerance 2 \
    --k 10 \
    --log_file ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_CellBC.mapping.log_middle10.txt

python 03_PAIR_Whitelist_Mapping.py \
    --input ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_CellBC.mapped_middle10.txt \
    --whitelist ${DATA_DIR}/code/PAIR_extraction/PAIR_perturb_whitelist_middle10.csv \
    --output ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped_middle10.txt \
    --k 15 \
    --tolerance 2 \
    --log ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped.log_middle10.txt

python 04_Generate_PAIR_UMI_matrix.py \
    -i ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped_middle10.txt \
    -s ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_sparse_UMI_matrix_middle10 \
    -f ${DATA_DIR}/PAIR_output_middle10/PAIR_matched_frequency_middle10.txt