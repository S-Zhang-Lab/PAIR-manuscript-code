import pandas as pd
import argparse
import logging
from time import time

# Set up argument parser
parser = argparse.ArgumentParser(description="Map rows from one file to a whitelist and extract matched rows.")
parser.add_argument("--input_file", required=True, help="Path to the input file (e.g., DDR1_unique_combinations_sorted.txt).")
parser.add_argument("--whitelist_file", required=True, help="Path to the whitelist file (e.g., Final_oligos_withID.csv).")
parser.add_argument("--output_file", required=True, help="Path to save the output file (e.g., DDR1_unique_combinations_sorted_mapped.txt).")

# Parse arguments
args = parser.parse_args()

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("mapping_log.log"),
        logging.StreamHandler()
    ]
)

def main():
    start_time = time()

    logging.info("Loading input file...")
    input_df = pd.read_csv(args.input_file, sep="\t")
    logging.info(f"Input file loaded. Total rows: {len(input_df)}")

    logging.info("Loading whitelist file...")
    whitelist_df = pd.read_csv(args.whitelist_file)
    logging.info(f"Whitelist file loaded. Total rows: {len(whitelist_df)}")

    logging.info("Performing mapping...")
    # Merge on the specified columns to ensure all three columns match
    merged_df = input_df.merge(
        whitelist_df,
        left_on=["Position1", "Position2_RC_RC", "BC14nt_RC"],
        right_on=["CRISPRa", "CasRx", "BC14"],
        how="inner"
    )
    logging.info(f"Mapping complete. Matched rows: {len(merged_df)}")

    logging.info("Rearranging columns for output...")
    # Arrange columns in the specified order
    output_columns = ["ID", "CRISPRa_name", "CasRx_name"] + list(input_df.columns)
    output_df = merged_df[output_columns]

    logging.info("Saving matched rows with rearranged columns...")
    output_df.to_csv(args.output_file, sep="\t", index=False)
    logging.info(f"Matched rows saved to {args.output_file}")

    total_rows = len(input_df)
    mapped_rows = len(output_df)
    percentage_mapped = (mapped_rows / total_rows) * 100

    logging.info(f"Total rows in input file: {total_rows}")
    logging.info(f"Rows mapped to whitelist: {mapped_rows}")
    logging.info(f"Percentage of rows mapped: {percentage_mapped:.2f}%")

    elapsed_time = time() - start_time
    logging.info(f"Processing completed in {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
