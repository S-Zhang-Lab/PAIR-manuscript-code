#!/bin/bash
#$ -m abe
#$ -pe smp 2 # reserving 2 cores
#$ -q long
#$ -N ref_genome


# Define reference path
REF_DIR="${PROJECT_ROOT:-.}/ref_genome"
cd "$REF_DIR" || { echo "Directory not found: $REF_DIR"; exit 1; }

# Download the latest human reference genome (≈11 GB)
wget -c https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz

# Extract the tarball
tar -xzvf refdata-gex-GRCh38-2024-A.tar.gz

# Optional: remove tar.gz to save space
# rm refdata-gex-GRCh38-2024-A.tar.gz

echo "✅ Download and extraction complete at $(date)"