#!/bin/bash
#$ -M dgan@nd.edu
#$ -m abe
#$ -pe smp 2  # Reserve 2 cores
#$ -q long
#$ -N ID1_TT  # job name
#$ -cwd  # run from current working directory

module load rclone

# Perform rclone copy
rclone copy "Dailin_box:UTSW32_XY_Naive_Omentum_XJ_Lib/01.RawData/UTSW32_ID1_SI_TT_G9" ./UTSW32_ID1_SI_TT_G9