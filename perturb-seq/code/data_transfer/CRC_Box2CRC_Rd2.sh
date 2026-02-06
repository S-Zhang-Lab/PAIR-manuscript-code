#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe
#$ -pe smp 2  # Reserve 2 cores
#$ -q long
#$ -N data_transfer  # job name

module load rclone

# Perform rclone copy
rclone copy "Dailin_box:UND_collabration/UTSW38_CC_perturb2/RAW" \
"/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/PAIRseq_Rd2" \
-P --create-empty-src-dirs --transfers=8 --checkers=8 \
--max-backlog 9999 --log-file=rclone_copy.log --log-level=INFO