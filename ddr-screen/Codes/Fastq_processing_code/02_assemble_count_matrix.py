import pandas as pd
import os
import argparse
from time import time

# Set up argument parser
parser = argparse.ArgumentParser(description="Assemble a count matrix from DDR files and whitelist.")
parser.add_argument("--whitelist_file", required=True, help="Path to the whitelist file (e.g., Final_oligos_withID.csv).")
parser.add_argument("--output_file", required=True, help="Path to save the assembled count matrix (e.g., count_matrix.csv).")
parser.add_argument("--input_files", nargs="+", required=True, help="List of DDR files to process (e.g., DDR1, DDR2, ...).")

# Parse arguments
args = parser.parse_args()

def assemble_count_matrix(whitelist_file, ddr_files, output_file):
    start_time = time()

    print("Loading whitelist file...")
    whitelist_df = pd.read_csv(whitelist_file)

    # Extract ID, CRISPRa, and CasRx columns
    count_matrix = whitelist_df[['ID', 'CRISPRa', 'CasRx', 'CRISPRa_name', 'CasRx_name']].copy()  # Explicitly create a copy
    count_matrix.set_index('ID', inplace=True)
    print(f"Whitelist file loaded. Total IDs: {len(count_matrix)}")

    for ddr_file in ddr_files:
        print(f"Processing {ddr_file}...")
        # Extract the DDR file name (e.g., DDR1, DDR2)
        column_name = os.path.basename(ddr_file).split("_")[0]

        # Load DDR file with comma as the delimiter
        ddr_df = pd.read_csv(ddr_file, sep=",")

        # Map counts to whitelist IDs
        ddr_counts = ddr_df.set_index('ID')['Count']
        count_matrix[column_name] = count_matrix.index.map(ddr_counts).fillna(0).astype(int)  # Use .loc for assignment
        print(f"Finished processing {ddr_file}.")

    print("Saving the count matrix...")
    count_matrix.to_csv(output_file)
    print(f"Count matrix saved to {output_file}")

    elapsed_time = time() - start_time
    print(f"Processing completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    assemble_count_matrix(args.whitelist_file, args.input_files, args.output_file)
