#!/bin/bash

# --- Configuration ---
# 1. Set your FASTQ file name
FASTQ_FILE="mRNA_TT_A12-SI_TT_A12_2333FLLT3_S2_L007_R1_001.fastq.gz"

# 2. Set the exact sequence you are searching for
TARGET_SEQUENCE="GTTTTAGAGCTAGGCCAACATGAGGATCACCCATGTCTGCAGGGCCTAGCAAGTTAAAATAAGGC"

# ---------------------

echo "Processing file: $FASTQ_FILE"
echo "Searching for sequence: $TARGET_SEQUENCE"
echo "..."

# 3. Run the single-pass awk command
zcat "$FASTQ_FILE" | awk -v seq="$TARGET_SEQUENCE" '
    BEGIN {
        total_reads = 0
        desired_reads = 0
    }
    
    # This selector runs only on sequence lines (line 2, 6, 10, etc.)
    (NR % 4 == 2) {
        # Count every sequence line as a total read
        total_reads++
        
        # Check if the sequence line ( $0 ) contains our target sequence ( seq )
        # index() checks for a literal string match.
        if (index($0, seq)) {
            # If it matches, increment the desired read counter
            desired_reads++
        }
    }
    
    # After processing the whole file, run this block
    END {
        if (total_reads > 0) {
            # Calculate the percentage
            percentage = (desired_reads / total_reads) * 100
            
            # Print the final report
            printf "Total Reads:     %d\n", total_reads
            printf "Desired Reads:   %d\n", desired_reads
            printf "Percentage:      %.4f%%\n", percentage
        } else {
            print "No reads found in file."
        }
    }
'