import gzip
from collections import defaultdict

# Define your CMO barcodes
barcodes = {
    "C0254_cmo": "AGTAAGTTCAGCGTA",
    "C0255_cmo": "AAGTATCGTTTCGCA",
    "C0256_cmo": "GGTTGCCAGATGTCA",
    "C0257_cmo": "TGTCTTTCCTGCCAG",
    "C0258_cmo": "CTCCTCTGCAATTAC",
    "C0259_cmo": "CAGTAGTCACGGTCA"
}

# FASTQ paths
R1_path = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/HTO_TN_A3/HTO_TN_A3-SI_TN_A3_2333FLLT3_S1_L007_R1_001.fastq.gz"
R2_path = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/HTO_TN_A3/HTO_TN_A3-SI_TN_A3_2333FLLT3_S1_L007_R2_001.fastq.gz"
output_path = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/cmo_counts_R1_R2.txt"

# Reverse complement helper
def revcomp(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

# Precompute reverse complements
rev_barcodes = {cmo: revcomp(seq) for cmo, seq in barcodes.items()}

# Initialize counters
counts_R1 = defaultdict(int)
counts_R2 = defaultdict(int)

def count_matches(fastq_path, seq_dict, counter):
    """Count number of reads containing each barcode in a gzipped FASTQ file."""
    with gzip.open(fastq_path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                read = line.strip()
                for cmo, barcode in seq_dict.items():
                    if barcode in read:
                        counter[cmo] += 1

# Count in R2 (forward sequences)
print("Counting in R2 (forward CMO sequences)...")
count_matches(R2_path, barcodes, counts_R2)

# Count in R1 (reverse complements)
print("Counting in R1 (reverse complement CMO sequences)...")
count_matches(R1_path, rev_barcodes, counts_R1)

# Save results
with open(output_path, "w") as out:
    out.write("CMO_ID\tForward_seq\tReverse_complement\tR2_count\tR1_count\n")
    for cmo in barcodes:
        out.write(f"{cmo}\t{barcodes[cmo]}\t{rev_barcodes[cmo]}\t{counts_R2[cmo]}\t{counts_R1[cmo]}\n")

print(f"✅ Done! Results saved to {output_path}")