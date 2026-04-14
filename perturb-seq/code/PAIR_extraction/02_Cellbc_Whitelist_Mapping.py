import argparse
import pandas as pd
from Levenshtein import distance
from collections import defaultdict

def log_to_file(log_file, message):
    """Log a message to both the console and a log file."""
    print(message)
    with open(log_file, 'a') as f:
        f.write(message + '\n')

def load_whitelist(file_path, log_file):
    """Load and reshape the CellBC whitelist, removing the "-1" suffix."""
    log_to_file(log_file, f"Loading whitelist from {file_path}...")
    whitelist_df = pd.read_csv(file_path, header=None)
    whitelist = whitelist_df[1].str.replace("-1", "", regex=False).tolist()
    log_to_file(log_file, f"Whitelist loaded: {len(whitelist)} entries.")
    log_to_file(log_file, "First few whitelist entries:")
    log_to_file(log_file, str(whitelist[:5]))
    return whitelist

def build_kmer_index(whitelist, k, log_file):
    """Build a k-mer index from the whitelist for quick lookup."""
    kmer_index = defaultdict(set)
    for ref_bc in whitelist:
        for i in range(len(ref_bc) - k + 1):
            kmer = ref_bc[i:i + k]
            kmer_index[kmer].add(ref_bc)
    log_to_file(log_file, f"K-mer index built with k-mer size {k}.")
    return kmer_index

def find_candidates(cell_bc, kmer_index, k):
    """Find candidate whitelist entries for a given CellBC using k-mers."""
    candidates = set()
    for i in range(len(cell_bc) - k + 1):
        kmer = cell_bc[i:i + k]
        if kmer in kmer_index:
            candidates.update(kmer_index[kmer])
    return candidates

def preprocess_for_distance(seq1, seq2):
    """Align sequences by ignoring positions with 'N' in seq1."""
    filtered_seq1 = ''.join([s1 for s1, s2 in zip(seq1, seq2) if s1 != 'N'])
    filtered_seq2 = ''.join([s2 for s1, s2 in zip(seq1, seq2) if s1 != 'N'])
    return filtered_seq1, filtered_seq2

def map_cellbc(cell_bc, kmer_index, whitelist, tolerance, k):
    """Map a single CellBC to the whitelist, allowing error correction."""
    candidates = find_candidates(cell_bc, kmer_index, k)
    for ref_bc in candidates:
        preprocessed_cell_bc, preprocessed_ref_bc = preprocess_for_distance(cell_bc, ref_bc)
        if distance(preprocessed_cell_bc, preprocessed_ref_bc) <= tolerance:
            return ref_bc
    return None

def map_cellbc_to_whitelist(extracted_data_path, whitelist, output_path, unmapped_path, tolerance, k, log_file):
    """Map CellBC values to the whitelist, allowing error correction."""
    log_to_file(log_file, f"Loading extracted data from {extracted_data_path}...")
    data = pd.read_csv(extracted_data_path, sep="\t", chunksize=1000)
    log_to_file(log_file, f"Processing data in chunks of 1000 rows.")

    log_to_file(log_file, f"Parameters: Tolerance = {tolerance}, K-mer size = {k}.")

    log_to_file(log_file, "Building k-mer index...")
    kmer_index = build_kmer_index(whitelist, k, log_file)

    log_to_file(log_file, "Starting CellBC mapping...")
    mapped_count = 0
    total_rows = 0

    with open(output_path, 'w') as mapped_file, open(unmapped_path, 'w') as unmapped_file:
        # Write headers to both files
        mapped_file.write("CellBC\tUMI\tCRISPRa_gRNA\tCasRx_crRNA\tMapped_CellBC\n")
        unmapped_file.write("CellBC\tUMI\tCRISPRa_gRNA\tCasRx_crRNA\n")

        for chunk in data:
            mapped_rows = []
            unmapped_rows = []

            for _, row in chunk.iterrows():
                total_rows += 1
                cell_bc = row["CellBC"]
                mapped_cellbc = map_cellbc(cell_bc, kmer_index, whitelist, tolerance, k)
                if mapped_cellbc:
                    mapped_count += 1
                    row["Mapped_CellBC"] = mapped_cellbc
                    mapped_rows.append(row)
                else:
                    unmapped_rows.append(row)

            # Write incrementally
            if mapped_rows:
                pd.DataFrame(mapped_rows).to_csv(mapped_file, sep="\t", header=False, index=False, mode='a')
            if unmapped_rows:
                pd.DataFrame(unmapped_rows).to_csv(unmapped_file, sep="\t", header=False, index=False, mode='a')

            log_to_file(log_file, f"Processed {total_rows} rows so far. {mapped_count} rows mapped.")

    # Final statistics
    mapped_percentage = (mapped_count / total_rows) * 100
    log_to_file(log_file, f"Mapping complete. {mapped_count} out of {total_rows} rows mapped ({mapped_percentage:.2f}%).")

def main():
    parser = argparse.ArgumentParser(description="Map CellBC values to a whitelist with error correction.")
    parser.add_argument("--whitelist", required=True, help="Path to the CellBC whitelist file.")
    parser.add_argument("--extracted_data", required=True, help="Path to the extracted data table.")
    parser.add_argument("--output", required=True, help="Path to the output file for mapped rows.")
    parser.add_argument("--unmapped", required=True, help="Path to the output file for unmapped rows.")
    parser.add_argument("--tolerance", type=int, default=1, help="Number of bit tolerance for error correction.")
    parser.add_argument("--k", type=int, default=6, help="K-mer size for fast mapping.")
    parser.add_argument("--log_file", required=True, help="Path to the log file.")

    args = parser.parse_args()

    # Log file setup
    with open(args.log_file, 'w') as f:
        f.write("Log for CellBC Mapping\n")

    # Load whitelist
    whitelist = load_whitelist(args.whitelist, args.log_file)

    # Map CellBC to whitelist
    map_cellbc_to_whitelist(args.extracted_data, whitelist, args.output, args.unmapped, args.tolerance, args.k, args.log_file)

if __name__ == "__main__":
    main()
