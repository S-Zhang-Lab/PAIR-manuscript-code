import gzip

def count_fastq_reads(fastq_path):
    """Count total number of reads in a FASTQ or FASTQ.GZ file."""
    open_func = gzip.open if fastq_path.endswith(".gz") else open
    with open_func(fastq_path, "rt") as f:
        total_lines = sum(1 for _ in f)
    return total_lines // 4  # 4 lines per read

# FASTQ paths
R1_path = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/HTO_TN_A3/HTO_TN_A3-SI_TN_A3_2333FLLT3_S1_L007_R1_001.fastq.gz"
R2_path = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/raw_data/01.RawData/HTO_TN_A3/HTO_TN_A3-SI_TN_A3_2333FLLT3_S1_L007_R2_001.fastq.gz"


R1_reads = count_fastq_reads(R1_path)
R2_reads = count_fastq_reads(R2_path)

print(f"R1 total reads: {R1_reads}")
print(f"R2 total reads: {R2_reads}")

# Optional: consistency check
if R1_reads != R2_reads:
    print("⚠️ Warning: R1 and R2 have different numbers of reads!")
else:
    print(f"✅ Total paired reads: {R1_reads}")