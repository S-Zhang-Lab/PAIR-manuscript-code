# PAIR-DDR

Investigation of Transcriptome Interaction in DNA Repair Pathway by Programmable CRISPR Paired Sequencing

### Background introduction (CC add)

### Design of Dual crRNA cassette

The dual regulatory crRNA cassette is with the following structure:

GAGGGCCTATTTCCCATGATTcgtctcacaccNNNNNNNNNNNNNNNNNNNNgttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAACnnnnnnnnnnnnnnnnnnnnnnnTTTTTTT[BC12]ctacagagacgcacttgtacttcagcggtca

N: 20 bases gRNA for CRISPRa

n: 23 bases crRNA for Cas13d

[BC12] : 12 digit BC with hamming distance \> 3, which allow one bit error correction
