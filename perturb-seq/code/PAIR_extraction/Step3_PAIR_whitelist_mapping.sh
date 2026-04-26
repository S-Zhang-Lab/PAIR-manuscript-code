#!/bin/bash
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_whitelist   # Specify job name

module load python
source ${PROJECT_ROOT:-.}/venv/bin/activate

DATA_DIR=${PROJECT_ROOT:-.}

python 03_PAIR_Whitelist_Mapping.py \
    --input ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.mapped.txt \
    --whitelist ${DATA_DIR}/code/PAIR_extraction/PAIR_perturb_whitelist.csv \
    --output ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped.txt \
    --k 15 \
    --tolerance 2 \
    --log ${DATA_DIR}/PAIR_output/PAIR_matched_extracted_data_table_CellBC.PAIR.mapped.log.txt
