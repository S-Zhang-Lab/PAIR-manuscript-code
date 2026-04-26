#!/bin/bash
#$ -m abe               # Send mail when job begins, ends and aborts
#$ -q long              # Specify queue
#$ -pe smp 4            # Specify number of cores to use.
#$ -N PAIR_redo       # Specify job name

module load python
source ${PROJECT_ROOT:-.}/venv/bin/activate

DATA_DIR=${PROJECT_ROOT:-.}/raw_data/01.RawData/PAIR_TT_C12
OUT_DIR=${PROJECT_ROOT:-.}/PAIR_output

# full run
python 01_Extract_PAIR.py   --input_r1 ${DATA_DIR}/PAIR_TT_C12-SI_TT_C12_2333FLLT3_S3_L007_R1_001.fastq.gz \
                            --input_r2 ${DATA_DIR}/PAIR_TT_C12-SI_TT_C12_2333FLLT3_S3_L007_R2_001.fastq.gz \
                            --output_r1 ${OUT_DIR}/PAIR_matched_R1_mismatch01_AdjR2.fastq.gz \
                            --output_r2 ${OUT_DIR}/PAIR_matched_R2_mismatch01_AdjR2.fastq.gz \
                            --output_table ${OUT_DIR}/PAIR_matched_extracted_data_table_mismatch01_AdjR2.txt \
                            --whitelist ${PROJECT_ROOT:-.}/code/PAIR_extraction/PAIR_perturb_whitelist.csv



