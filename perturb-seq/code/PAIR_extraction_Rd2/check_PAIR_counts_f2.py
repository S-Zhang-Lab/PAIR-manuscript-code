import gzip
import csv

# --------------------------
# Configurable file paths
# --------------------------
R1_path = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/PAIRseq_Rd2/01.RawData/CC_PAIR_f2/CC_PAIR_f2-SI_TT_F2_23357TLT4_S6_L004_R1_001.fastq.gz"
R2_path = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/PAIRseq_Rd2/01.RawData/CC_PAIR_f2/CC_PAIR_f2-SI_TT_F2_23357TLT4_S6_L004_R2_001.fastq.gz"
guide_csv = "/afs/crc.nd.edu/group/StatDataMine/dm008/Dailin_Gan/Siyuan/PAIR-perturb-seq_2025-10/code/PAIR_extraction_Rd2/PAIR_perturb_whitelist.csv"

# --------------------------
# Helper functions
# --------------------------
def open_fastq(path):
    """Open gzip or normal FASTQ in text mode."""
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    complement = str.maketrans("ATCGNatcgn", "TAGCNtagcn")
    return seq.translate(complement)[::-1]

def read_fastq_sequences(path):
    """Yield every 2nd line (sequence line) from FASTQ file."""
    with open_fastq(path) as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()      # sequence
            f.readline()                    # '+'
            f.readline()                    # quality
            yield seq

# --------------------------
# Load guide sequences
# --------------------------
crispr_a_guides = {}
casrx_guides = {}

with open(guide_csv) as f:
    reader = csv.DictReader(f)
    for row in reader:
        crispr_a_guides[row["CRISPRa"]] = row["CRISPRa_name"]
        casrx_guides[row["CasRx"]] = row["CasRx_name"]

# Initialize counters
crispr_a_counts = {seq: 0 for seq in crispr_a_guides}
casrx_counts = {seq: 0 for seq in casrx_guides}
total_R1_reads = 0
total_R2_reads = 0

# --------------------------
# Count CRISPRa matches in R1
# --------------------------
for seq in read_fastq_sequences(R1_path):
    total_R1_reads += 1
    for guide_seq in crispr_a_counts:
        if guide_seq in seq:  # substring match
            crispr_a_counts[guide_seq] += 1

# --------------------------
# Count CasRx matches in R2 (reverse complement)
# --------------------------
for seq in read_fastq_sequences(R2_path):
    total_R2_reads += 1
    rc_seq = reverse_complement(seq)
    for guide_seq in casrx_counts:
        if guide_seq in rc_seq:  # substring match in reverse complement
            casrx_counts[guide_seq] += 1

# --------------------------
# Reporting
# --------------------------
print("\n=== CRISPRa (R1 substring matches) ===")
for seq, count in crispr_a_counts.items():
    print(f"{crispr_a_guides[seq]:<15} {seq:<25} count = {count}")
print(f"Total CRISPRa matches: {sum(crispr_a_counts.values())}")
print(f"Total R1 reads: {total_R1_reads}")
print(f"Match rate: {100 * sum(crispr_a_counts.values()) / total_R1_reads:.5f}%")

print("\n=== CasRx (R2 reverse complement substring matches) ===")
for seq, count in casrx_counts.items():
    print(f"{casrx_guides[seq]:<15} {seq:<25} count = {count}")
print(f"Total CasRx matches: {sum(casrx_counts.values())}")
print(f"Total R2 reads: {total_R2_reads}")
print(f"Match rate: {100 * sum(casrx_counts.values()) / total_R2_reads:.5f}%")

if total_R1_reads == total_R2_reads:
    print(f"\n✅ Total paired reads: {total_R1_reads}")
else:
    print(f"\n⚠️ Warning: R1 ({total_R1_reads}) and R2 ({total_R2_reads}) read counts differ.")