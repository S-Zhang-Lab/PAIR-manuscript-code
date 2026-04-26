#!/bin/bash
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_CellBC   # Specify job name

module load python
source ${PROJECT_ROOT:-.}/venv/bin/activate

DATA_DIR=${PROJECT_ROOT:-.}

python 02_Cellbc_Whitelist_Mapping.py \
    --whitelist ${DATA_DIR}/code/PAIR_extraction/mRNA_barcodes_with_prefix.csv \
    --extracted_data ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_mismatch01_AdjR2.txt \
    --output ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.mapped.txt \
    --unmapped ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.unmapped.txt \
    --tolerance 2 \
    --k 10 \
    --log_file ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.mapping.log.txt

