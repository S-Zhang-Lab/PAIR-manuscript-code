import pandas as pd
import logging
import argparse
from typing import List, Dict, Any

def hamming_distance(str1: str, str2: str) -> int:
    """Calculates the Hamming distance between two strings."""
    distance = 0
    for c1, c2 in zip(str1, str2):
        if c1 != 'N' and c2 != 'N' and c1 != c2:
            distance += 1
    return distance

def map_to_whitelist(input_df: pd.DataFrame, whitelist_df: pd.DataFrame, tolerance: int) -> pd.DataFrame:
    """Maps sequences from input_df to whitelist_df based on k-mer matching and Hamming distance."""
    logging.info("Starting k-mer mapping process...")
    total_sequences = len(input_df)
    mapped_sequences = 0

    k = 20  # k-mer size

    whitelist_index: Dict[str, List[int]] = {}
    for i, row in whitelist_df.iterrows():
        seq = row['CRISPRa'] + row['CasRx'] + row['BC14']
        for j in range(len(seq) - k + 1):
            kmer = seq[j:j + k]
            if kmer not in whitelist_index:
                whitelist_index[kmer] = []
            whitelist_index[kmer].append(i)

    mapped_data: List[Dict[str, Any]] = []
    for idx, row in input_df.iterrows():
        if idx % 5000 == 0:
            logging.info(f"Processed {idx}/{total_sequences} sequences...")

        query = row['Position1'] + row['Position2_RC_RC'] + row['BC14nt_RC']
        matched_whitelist_indices = set()

        for i in range(len(query) - k + 1):
            kmer = query[i:i + k]
            if kmer in whitelist_index:
                matched_whitelist_indices.update(whitelist_index[kmer])

        best_match = None
        min_distance = float('inf')

        for idx in matched_whitelist_indices:
            entry = whitelist_df.iloc[idx]
            whitelist_seq = entry['CRISPRa'] + entry['CasRx'] + entry['BC14']

            dist = hamming_distance(whitelist_seq, query)
            if dist <= tolerance and dist < min_distance:
                min_distance = dist
                best_match = entry

        if best_match is not None:
            mapped_sequences += 1
            mapped_data.append({
                'Position1': row['Position1'],
                'Position2_RC_RC': row['Position2_RC_RC'],
                'BC14nt_RC': row['BC14nt_RC'],
                'Count': row['Count'],
                'ID': best_match['ID'],
                'CRISPRa': best_match['CRISPRa'],
                'CRISPRa_name': best_match['CRISPRa_name'],
                'CasRx': best_match['CasRx'],
                'CasRx_name': best_match['CasRx_name'],
                'BC14': best_match['BC14'],
                'Hamming_Distance': min_distance
            })

    mapping_rate = (mapped_sequences / total_sequences) * 100
    logging.info(f"Mapping complete: {mapped_sequences}/{total_sequences} sequences mapped ({mapping_rate:.2f}%)")

    return pd.DataFrame(mapped_data)

def collapse_and_merge(mapped_df: pd.DataFrame, whitelist_df: pd.DataFrame) -> pd.DataFrame:
    """Collapses the mapped DataFrame, merges it with the whitelist, and sorts the result."""
    logging.info("Collapsing and merging data...")
    
    # Collapse mapped_df by 'ID', summing the 'Count' column
    df1_collapsed = mapped_df.groupby('ID')['Count'].sum().reset_index()

    # Find the first row with Hamming_Distance = 0 for each ID
    df1_unique = mapped_df[mapped_df['Hamming_Distance'] == 0].drop_duplicates(subset=['ID'], keep='first')

    # Merge df1_collapsed and whitelist_df on the 'ID' column
    merged_df = pd.merge(df1_collapsed, whitelist_df, on='ID', how='inner')

    # Select columns 'CRISPRa_name' and 'CasRx_name' from df1_unique
    df1_selected = df1_unique[['ID', 'CRISPRa_name', 'CasRx_name']]

    # Merge merged_df and df1_selected on the 'ID' column
    merged_df = pd.merge(merged_df, df1_selected, on='ID', how='left')

    # Sort the merged DataFrame by `Count` in descending order
    final_sorted_df = merged_df.sort_values(by='Count', ascending=False)

    return final_sorted_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Combined K-mer mapping and collapse script.")
    parser.add_argument('--input_file', type=str, required=True, help="Input tab-delimited file with sequences to map.")
    parser.add_argument('--whitelist_file', type=str, required=True, help="Whitelist CSV file with reference sequences.")
    parser.add_argument('--output_file', type=str, required=True, help="Output CSV file for final results.")
    parser.add_argument('--tolerance', type=int, required=True, help="Hamming distance tolerance.")
    parser.add_argument('--count_threshold', type=int, required=True, help="Minimum count threshold for filtering.")

    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    try:
        logging.info("Starting script execution...")
        # Read input and whitelist files
        logging.info(f"Reading input file: {args.input_file}")
        input_df = pd.read_csv(args.input_file, sep="\t")
        logging.info(f"Reading whitelist file: {args.whitelist_file}")
        whitelist_df = pd.read_csv(args.whitelist_file)

        # Filter input data by count threshold
        logging.info(f"Filtering input data with count threshold: {args.count_threshold}")
        input_df_filtered = input_df[input_df['Count'] >= args.count_threshold]
        logging.info(f"Filtered input data. Remaining sequences: {len(input_df_filtered)}")

        # Perform mapping
        mapped_df = map_to_whitelist(input_df_filtered, whitelist_df, args.tolerance)

        # Collapse, merge, and sort
        final_df = collapse_and_merge(mapped_df, whitelist_df)

        # Save final results
        logging.info(f"Saving final data to: {args.output_file}")
        final_df.to_csv(args.output_file, index=False)
        logging.info("Script execution completed successfully.")

    except Exception as e:
        logging.error(f"An error occurred: {str(e)}")
        raise