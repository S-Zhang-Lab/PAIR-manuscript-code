#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_Steps2to4_f2   # Specify job name

set -e  # stop on errors

module load python
source /afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/venv/bin/activate

DATA_DIR=/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10

python 02_Cellbc_Whitelist_Mapping.py \
    --whitelist ${DATA_DIR}/code/PAIR_extraction_Rd2/f2_combined_sample_filtered_barcodes.csv \
    --extracted_data ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_f2_table_mismatch01_AdjR2.txt \
    --output ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_extracted_data_table_CellBC.mapped.txt \
    --unmapped ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_extracted_data_table_CellBC.unmapped.txt \
    --tolerance 2 \
    --k 10 \
    --log_file ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_extracted_data_table_CellBC.mapping.log.txt

python 03_PAIR_Whitelist_Mapping.py \
    --input ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_extracted_data_table_CellBC.mapped.txt \
    --whitelist ${DATA_DIR}/code/PAIR_extraction_Rd2/PAIR_perturb_whitelist.csv \
    --output ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped.txt \
    --k 15 \
    --tolerance 2 \
    --log ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped.log.txt

python 04_Generate_PAIR_UMI_matrix.py \
    -i ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped.txt \
    -s ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_sparse_UMI_matrix \
    -f ${DATA_DIR}/PAIRseq_Rd2/PAIR_f2/PAIR_matched_frequency.txt