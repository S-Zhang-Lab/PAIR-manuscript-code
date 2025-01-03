#!/usr/bin/env python3

import gzip
import re
import argparse
import logging
import sys
from datetime import datetime

def setup_logging(log_file, verbose):
    """
    Sets up the logging configuration.
    Logs are written to both the console and a log file.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    
    # Formatter with timestamp
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
        description="Extract specific sequences from paired-end FASTQ files based on defined patterns."
    )
    parser.add_argument(
        "-r1", "--read1",
        required=True,
        help="Path to the R1 FASTQ.gz file."
    )
    parser.add_argument(
        "-r2", "--read2",
        required=True,
        help="Path to the R2 FASTQ.gz file."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output TXT file."
    )
    parser.add_argument(
        "-l", "--log",
        default=f"extract_sequences_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
        help="Path to the log file. Defaults to 'extract_sequences_YYYYMMDD_HHMMSS.log'."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging to console."
    )
    parser.add_argument(
        "-p", "--progress",
        type=int,
        default=100000,
        help="Number of read pairs to process before logging progress. Default is 100,000."
    )
    return parser.parse_args()

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    
    Args:
        seq (str): DNA sequence.
    
    Returns:
        str: Reverse complement of the input sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                  'N': 'N', 'n': 'n'}
    try:
        # Reverse the sequence and complement each base
        rc_seq = ''.join([complement.get(base, 'N') for base in reversed(seq)])
        return rc_seq.upper()
    except Exception as e:
        # If unexpected characters are present, log and return 'N's
        return 'N' * len(seq)

def extract_sequences_with_header_matching(r1_path, r2_path, output_path, logger, progress_interval=100000):
    """
    Extracts Position1 from R1 and BC14nt & Position2_RC from R2 based on specified patterns.
    Ensures that R1 and R2 reads are correctly paired by matching headers.
    Converts BC14nt and Position2_RC to their reverse complements.
    Writes the extracted data to a text file and logs QC statistics.
    """
    
    # Define patterns
    pattern_r1 = re.compile(
        r'GTGGAAAGGACGAAACACCG(?P<position1>.{20})GTTTTAGAGCTAGGCCAACATGAGGATCACCCATGTCTGCAG'
    )
    
    pattern_r2 = re.compile(
        r'TGACCGCTGAAGTACAAGTGGTAGAGTAG(?P<bc14nt>.{14})A{7}(?P<position2_rc>.{23})GTTTCAAACCCC'
    )
    
    total_reads = 0
    successful_reads = 0
    header_mismatches = 0
    
    try:
        with gzip.open(r1_path, 'rt') as r1_file, gzip.open(r2_path, 'rt') as r2_file, open(output_path, 'w') as out_file:
            # Write header to output file
            out_file.write("Position1\tPosition2_RC_RC\tBC14nt_RC\n")
            logger.info(f"Started processing:\nR1: {r1_path}\nR2: {r2_path}\nOutput: {output_path}")
            
            while True:
                # Read four lines for R1
                r1_header = r1_file.readline()
                if not r1_header:
                    break  # End of file
                r1_seq = r1_file.readline().strip()
                r1_plus = r1_file.readline()
                r1_qual = r1_file.readline()
                
                # Read four lines for R2
                r2_header = r2_file.readline()
                if not r2_header:
                    logger.warning("R2 file ended before R1 file.")
                    break  # End of file
                r2_seq = r2_file.readline().strip()
                r2_plus = r2_file.readline()
                r2_qual = r2_file.readline()
                
                total_reads += 1
                
                # Extract read identifiers (everything before the first space)
                r1_id = r1_header.split()[0]
                r2_id = r2_header.split()[0]
                
                if r1_id != r2_id:
                    header_mismatches += 1
                    logger.warning(f"Header mismatch at read pair {total_reads}:")
                    logger.warning(f"R1 Header: {r1_header.strip()}")
                    logger.warning(f"R2 Header: {r2_header.strip()}")
                    continue  # Skip this pair
                
                # Search for patterns
                match_r1 = pattern_r1.search(r1_seq)
                match_r2 = pattern_r2.search(r2_seq)
                
                if match_r1 and match_r2:
                    position1 = match_r1.group('position1')
                    bc14nt = match_r2.group('bc14nt')
                    position2_rc = match_r2.group('position2_rc')
                    
                    # Convert BC14nt and Position2_RC to their reverse complements
                    bc14nt_rc = reverse_complement(bc14nt)
                    position2_rc_rc = reverse_complement(position2_rc)
                    
                    # Write to output file with reverse complemented sequences
                    out_file.write(f"{position1}\t{position2_rc_rc}\t{bc14nt_rc}\n")
                    
                    successful_reads += 1
                
                # Log progress at specified intervals
                if total_reads % progress_interval == 0:
                    logger.info(f"Processed {total_reads} read pairs. Successful extractions: {successful_reads}.")
    
    except Exception as e:
        logger.error("An error occurred during processing.")
        logger.exception(e)
        sys.exit(1)
    
    # Calculate percentage
    if total_reads > 0:
        success_percentage = (successful_reads / total_reads) * 100
        mismatch_percentage = (header_mismatches / total_reads) * 100
    else:
        success_percentage = 0.0
        mismatch_percentage = 0.0
    
    # Final QC statistics
    logger.info("Processing completed.")
    logger.info("QC Statistics:")
    logger.info(f"Total read pairs processed: {total_reads}")
    logger.info(f"Successfully extracted reads: {successful_reads} ({success_percentage:.2f}%)")
    if header_mismatches > 0:
        logger.info(f"Header mismatches detected: {header_mismatches} ({mismatch_percentage:.2f}%)")
    else:
        logger.info("No header mismatches detected.")

def main():
    args = parse_arguments()
    
    # Setup logging
    logger = setup_logging(args.log, args.verbose)
    
    # Log the start of the script
    logger.info("Script started.")
    logger.info(f"Input R1 file: {args.read1}")
    logger.info(f"Input R2 file: {args.read2}")
    logger.info(f"Output file: {args.output}")
    logger.info(f"Log file: {args.log}")
    
    # Start extraction
    extract_sequences_with_header_matching(
        r1_path=args.read1,
        r2_path=args.read2,
        output_path=args.output,
        logger=logger,
        progress_interval=args.progress
    )
    
    logger.info("Script finished successfully.")

if __name__ == "__main__":
    main()
