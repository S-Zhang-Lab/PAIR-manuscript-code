import pandas as pd
import argparse
import logging
from scipy.sparse import csr_matrix, save_npz
from scipy.io import mmwrite

def process_input_file(input_file, sparse_output_base, freq_output_file):
    """
    Process the input file to calculate frequency of Mapped_name,
    generate a sparse UMI matrix, and create output files.

    Parameters:
        input_file (str): Path to the input file.
        sparse_output_base (str): Base name for the sparse matrix output files (no extension).
        freq_output_file (str): Path to the output frequency table file.
    """
    # Initialize logging
    logging.info(f"Starting processing of {input_file}")
    
    # Step 1: Load the data
    try:
        data = pd.read_csv(input_file, sep="\t")
        logging.info(f"Input file {input_file} loaded successfully with {len(data)} rows.")
    except Exception as e:
        logging.error(f"Failed to load input file {input_file}: {e}")
        return

    # Step 2: Directly use "Mapped_name" as the unique gRNA combination
    data['gRNA_combination'] = data['Mapped_name']
    logging.info(f"Used 'Mapped_name' as the unique gRNA combination.")

    # Step 3: Calculate frequency of "Mapped_name"
    try:
        frequency = data['Mapped_name'].value_counts().reset_index()
        frequency.columns = ['Mapped_name', 'Frequency']
        frequency_sorted = frequency.sort_values(by='Frequency', ascending=False)
        frequency_sorted.to_csv(freq_output_file, sep="\t", index=False)
        logging.info(f"Frequency table saved to {freq_output_file} with {len(frequency_sorted)} rows.")
    except Exception as e:
        logging.error(f"Failed to calculate and save frequency table: {e}")
        return

    # Step 4: Group by Mapped_CellBC and gRNA_combination and count unique UMIs
    try:
        grouped = (
            data.groupby(['Mapped_CellBC', 'gRNA_combination'])['UMI']
            .nunique()
            .reset_index()
            .rename(columns={'UMI': 'UMI_count'})
        )
        logging.info(f"Grouped data into {grouped.shape[0]} unique (Mapped_CellBC, gRNA_combination) pairs.")
    except Exception as e:
        logging.error(f"Failed to group data for UMI matrix generation: {e}")
        return

    # Step 5: Pivot the data to create a count matrix
    try:
        count_matrix = grouped.pivot(index='Mapped_CellBC', columns='gRNA_combination', values='UMI_count')
        count_matrix.fillna(0, inplace=True)
        logging.info("Pivoted data into a matrix format.")
    except Exception as e:
        logging.error(f"Failed to pivot data into matrix format: {e}")
        return

    # Step 6: Convert to sparse matrix and save in MTX and NPZ formats
    try:
        sparse_matrix = csr_matrix(count_matrix.values)

        # Save the sparse matrix in Matrix Market format
        mtx_output_file = sparse_output_base + ".mtx"
        mmwrite(mtx_output_file, sparse_matrix)
        logging.info(f"Matrix Market file saved to {mtx_output_file}")

        # Save row and column labels for reference
        count_matrix.index.to_series().to_csv(sparse_output_base + "_rows.txt", index=False, header=False)
        count_matrix.columns.to_series().to_csv(sparse_output_base + "_columns.txt", index=False, header=False)
        logging.info(f"Row labels saved to {sparse_output_base}_rows.txt")
        logging.info(f"Column labels saved to {sparse_output_base}_columns.txt")

        # Save matrix in .npz format for Python compatibility
        npz_output_file = sparse_output_base + ".npz"
        save_npz(npz_output_file, sparse_matrix)
        logging.info(f"Sparse matrix also saved to {npz_output_file}")

        # Log statistics
        logging.info(f"Matrix statistics: {sparse_matrix.shape[0]} cells, {sparse_matrix.shape[1]} gRNA combinations.")
        logging.info(f"Non-zero entries in matrix: {sparse_matrix.nnz}")
    except Exception as e:
        logging.error(f"Failed to save sparse matrix files: {e}")
        return


if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Process PAIR data and calculate Mapped_name frequencies.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input file.")
    parser.add_argument("-s", "--sparse_output", required=True, help="Base name for the sparse matrix output files (no extension).")
    parser.add_argument("-f", "--freq_output", required=True, help="Path to the output frequency table file.")
    args = parser.parse_args()

    # Set up logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler("PAIR_matched_sparse_UMI_matrix.log"),
            logging.StreamHandler()
        ]
    )

    # Run the processing function
    process_input_file(args.input, args.sparse_output, args.freq_output)
