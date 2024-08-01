# PAIR-DDR

Investigation of Transcriptome Interaction in DNA Repair Pathway by Programmable CRISPR Paired Sequencing

### Background introduction (CC add)

### Design of Dual crRNA cassette

The dual regulatory crRNA cassette is with the following structure:

<<<<<<< HEAD
GAGGGCCTATTTCCCATGATTcgtctcacaccNNNNNNNNNNNNNNNNNNNNgttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAACnnnnnnnnnnnnnnnnnnnnnnnTTTTTTT[BC10]ctacagagacgcacttgtacttcagcggtca
=======
GAGGGCCTATTTCCCATGATTcgtctcacaccgNNNNNNNNNNNNNNNNNNNNgttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAACnnnnnnnnnnnnnnnnnnnnnnnTTTTTTTctacagagacgcacttgtacttcagcggtca
>>>>>>> 1b8b3f954286b3757eb894e6303ed57be6d26fbf

N: 20 bases gRNA for CRISPRa
n: 23 bases crRNA for Cas13d
[BC10] : 10 digit BC with hamming distance \>4
