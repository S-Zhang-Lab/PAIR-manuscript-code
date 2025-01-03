#!/usr/bin/env python3

import pandas as pd
import argparse
import logging
import sys
from datetime import datetime

def setup_logging(log_file, verbose):
    """
    Sets up logging to console and a log file.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    # Formatter
    formatter = logging.Formatter(
        fmt='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # File handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    return logger

def parse_arguments():
    """
    Parses command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Calculate frequency of unique Position1 + Position2_RC_RC + BC14nt combinations."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the extracted sequences TXT file (e.g., DDR1_extraction_output.txt)."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output file to save sorted unique combinations."
    )
    parser.add_argument(
        "-l", "--log",
        default=f"calculate_frequencies_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
        help="Path to the log file. Defaults to 'calculate_frequencies_YYYYMMDD_HHMMSS.log'."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging to console."
    )
    return parser.parse_args()

def calculate_frequencies_pandas(input_file, output_file, logger):
    """
    Calculates frequencies using Pandas and saves the sorted result.
    """
    try:
        logger.info(f"Reading input file: {input_file}")
        # Read the input file
        df = pd.read_csv(input_file, sep='\t')
        
        logger.info("Calculating frequencies of unique combinations.")
        # Group by the three columns and count occurrences
        frequency_df = df.groupby(['Position1', 'Position2_RC_RC', 'BC14nt_RC']).size().reset_index(name='Count')
        
        logger.info("Sorting the combinations by frequency in descending order.")
        # Sort by Count descending
        frequency_df_sorted = frequency_df.sort_values(by='Count', ascending=False)
        
        logger.info(f"Saving the sorted frequencies to: {output_file}")
        # Save to output file
        frequency_df_sorted.to_csv(output_file, sep='\t', index=False)
        
        logger.info("Frequency calculation and sorting completed successfully.")
    
    except Exception as e:
        logger.error("An error occurred during frequency calculation.")
        logger.exception(e)
        sys.exit(1)

def calculate_frequencies_standard(input_file, output_file, logger):
    """
    Calculates frequencies using standard Python libraries and saves the sorted result.
    """
    from collections import defaultdict
    
    frequency_dict = defaultdict(int)
    
    try:
        logger.info(f"Reading input file: {input_file}")
        with open(input_file, 'r') as infile:
            header = infile.readline()  # Skip header
            for line_num, line in enumerate(infile, start=2):
                line = line.strip()
                if not line:
                    continue  # Skip empty lines
                parts = line.split('\t')
                if len(parts) != 3:
                    logger.warning(f"Line {line_num} is malformed: {line}")
                    continue
                key = tuple(parts)  # (Position1, Position2_RC_RC, BC14nt_RC)
                frequency_dict[key] += 1
                
                if line_num % 100000 == 0:
                    logger.info(f"Processed {line_num} lines.")
        
        logger.info("Calculating frequencies of unique combinations.")
        # Convert the frequency dictionary to a list of tuples
        frequency_list = [(*key, count) for key, count in frequency_dict.items()]
        
        logger.info("Sorting the combinations by frequency in descending order.")
        # Sort the list by count descending
        frequency_list_sorted = sorted(frequency_list, key=lambda x: x[3], reverse=True)
        
        logger.info(f"Saving the sorted frequencies to: {output_file}")
        # Write to the output file
        with open(output_file, 'w') as outfile:
            outfile.write("Position1\tPosition2_RC_RC\tBC14nt_RC\tCount\n")
            for item in frequency_list_sorted:
                outfile.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\n")
        
        logger.info("Frequency calculation and sorting completed successfully.")
    
    except Exception as e:
        logger.error("An error occurred during frequency calculation.")
        logger.exception(e)
        sys.exit(1)

def main():
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logging(args.log, args.verbose)
    
    # Log the start of the script
    logger.info("Script started.")
    logger.info(f"Input file: {args.input}")
    logger.info(f"Output file: {args.output}")
    
    # Choose which method to use
    # Uncomment the desired method
    
    # Method 1: Using Pandas
    calculate_frequencies_pandas(args.input, args.output, logger)
    
    # Method 2: Using Standard Python Libraries
    # calculate_frequencies_standard(args.input, args.output, logger)
    
    logger.info("Script finished successfully.")

if __name__ == "__main__":
    main()
