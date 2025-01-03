# PAIR-DDR

Investigation of Transcriptome Interaction in DNA Repair Pathway by Programmable CRISPR Paired Sequencing

### Background introduction (CC add)

U6-driven PAIR RNA cassette structure can be find here: <https://benchling.com/s/seq-SQ0Ur1pqXan8cwLW510g?m=slm-h2zfzbW6z9kpwVhFuXMc>

### Design of Dual crRNA cassette

The dual regulatory crRNA cassette is with the following structure:

e1 \<- "AGGGCCTATTTCCCATGATTcgtctcacaccg"

e2 \<- "N20": CRISPRa gRNA

e3 \<- "gttttagagctaggccaacatgaggatcacccatgtctgcagggcctagcaagttaaaataaggctagtccgttatcaacttggccaacatgaggatcacccatgtctgcagggccaagtggcaccgagtcggtgcttCAAGTAAACCCCTACCAACTGGTCGGGGTTTGAAAC"

e4 \<- "n23": Cas13d (CasRx) crRNA

e5 \<- "TTTTTTT"

e6 \<- "BC14": 14 digit BC with hamming distance more than 5, which allow 2 bit error correction

e7 \<- "ctacagagacgcacttgtacttcagcggtc"
