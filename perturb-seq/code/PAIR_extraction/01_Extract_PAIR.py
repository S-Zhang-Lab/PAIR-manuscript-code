import gzip
# import re
import regex
from Bio.Seq import Seq
import argparse
from datetime import datetime


def log_message(message, log_file):
    """
    Write a log message with a timestamp to the console and a log file.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_entry = f"[{timestamp}] {message}"
    print(log_entry)
    with open(log_file, 'a') as log:
        log.write(log_entry + '\n')


def extract_matched_reads(input_r1, input_r2, output_r1, output_r2, r1_pattern, r2_pattern, log_file):
    """
    Extract reads matching specific patterns from R1 and R2 files.
    Writes data in chunks to minimize I/O overhead while maintaining detailed logging.
    """
    total_reads = 0
    matched_reads = 0
    write_buffer_r1 = []
    write_buffer_r2 = []
    buffer_size = 100000  # Adjust buffer size based on available memory

    with gzip.open(input_r1, 'rt') as r1, gzip.open(input_r2, 'rt') as r2, \
         gzip.open(output_r1, 'wt') as out_r1, gzip.open(output_r2, 'wt') as out_r2:

        while True:
            # Read four lines (FASTQ format)
            r1_lines = [r1.readline().strip() for _ in range(4)]
            r2_lines = [r2.readline().strip() for _ in range(4)]

            if not r1_lines[0] or not r2_lines[0]:
                break

            total_reads += 1

            # Check for matching patterns in R1 and R2
            # Allowing up to 1 mismatch using regex module
            if (
                regex.match(f".{{28}}(?:{regex.escape(r1_pattern)}){{e<=1}}", r1_lines[1][:42])
                and regex.search(f"(?:{regex.escape(r2_pattern)}){{e<=1}}", r2_lines[1])
            ):
                matched_reads += 1
                write_buffer_r1.append('\n'.join(r1_lines) + '\n')
                write_buffer_r2.append('\n'.join(r2_lines) + '\n')

            # Write to disk when buffer is full
            if len(write_buffer_r1) >= buffer_size:
                out_r1.writelines(write_buffer_r1)
                out_r2.writelines(write_buffer_r2)
                write_buffer_r1.clear()
                write_buffer_r2.clear()

            # Log detailed progress every 100,000 reads
            if total_reads % 100000 == 0:
                log_message(f"Processed {total_reads} reads so far...", log_file)
                log_message(f"Matched reads so far: {matched_reads}", log_file)

        # Write any remaining data in buffers
        if write_buffer_r1:
            out_r1.writelines(write_buffer_r1)
            out_r2.writelines(write_buffer_r2)

    log_message(f"Total reads processed: {total_reads}", log_file)
    log_message(f"Matched reads: {matched_reads}", log_file)
    if total_reads > 0:
        log_message(f"Percentage of matched reads: {matched_reads / total_reads * 100:.2f}%", log_file)
    else:
        log_message("No reads processed.", log_file)

def process_reads(input_r1, input_r2, output_table, whitelist, log_file):
    """
    Process matched reads and extract CellBC, UMI, CRISPRa_gRNA, and CasRx_crRNA.
    Writes results incrementally to minimize memory usage while maintaining detailed logging.
    """
    r1_pattern = "TTTCTTATATGGGG"
    # r2_pattern = "TTAAAGCGTTTCAAACCCCGACCAGTTGGTAGGGGTTTACTTG"
    r2_pattern = "TTGCTAGGACCGGCCTTAAAGC" # use the CS1 region pattern to capture Cas13d crRNA sequences

    processed_reads = 0
    matched_records = 0  # Initialize matched_records

    with gzip.open(input_r1, 'rt') as r1, gzip.open(input_r2, 'rt') as r2, open(output_table, 'w') as out_table:
        # Write header to the output file
        out_table.write("CellBC\tUMI\tCRISPRa_gRNA\tCasRx_crRNA\n")

        while True:
            # Read four lines (FASTQ format)
            r1_lines = [r1.readline().strip() for _ in range(4)]
            r2_lines = [r2.readline().strip() for _ in range(4)]

            if not r1_lines[0] or not r2_lines[0]:
                break

            processed_reads += 1

            # Log progress every 100,000 reads
            if processed_reads % 100000 == 0:
                log_message(f"Processed {processed_reads} matched reads so far...", log_file)

            r1_seq = r1_lines[1]
            r2_seq = r2_lines[1]

            # Extract components
            cell_bc = r1_seq[:16]
            umi = r1_seq[16:28]
            crispr_g_rna_start = r1_seq.find(r1_pattern) + len(r1_pattern)
            crispr_g_rna = r1_seq[crispr_g_rna_start:crispr_g_rna_start + 20]
            casrx_crrna_start = r2_seq.find(r2_pattern) + len(r2_pattern) + 36 # 36 is the length of RfxCas13d DR36 region
            casrx_crrna = r2_seq[casrx_crrna_start:casrx_crrna_start + 23] # Have questions about the length 23 or 22? Ans: 23
            casrx_crrna_rc = str(Seq(casrx_crrna).reverse_complement())

            out_table.write(f"{cell_bc}\t{umi}\t{crispr_g_rna}\t{casrx_crrna_rc}\n")
            matched_records += 1

            # Flush every 100,000 matched records
            if matched_records % 100000 == 0:
                out_table.flush()
                log_message(f"Flushed {matched_records} matched records to the table so far...", log_file)

        log_message(f"Processing of matched reads complete. Total matched reads processed: {processed_reads}", log_file)

def main():
    parser = argparse.ArgumentParser(description="Process paired-end sequencing data.")
    parser.add_argument("--input_r1", required=True, help="Input R1 FASTQ file (gzipped).")
    parser.add_argument("--input_r2", required=True, help="Input R2 FASTQ file (gzipped).")
    parser.add_argument("--output_r1", required=True, help="Output R1 matched FASTQ file (gzipped).")
    parser.add_argument("--output_r2", required=True, help="Output R2 matched FASTQ file (gzipped).")
    parser.add_argument("--output_table", required=True, help="Output table with extracted data.")
    parser.add_argument("--whitelist", required=True, help="sgRNA whitelist file.")

    args = parser.parse_args()

    log_file = "01_Extract_PAIR_mismatch01.log"

    log_message("Starting extraction of matched reads...", log_file)
    extract_matched_reads(args.input_r1, args.input_r2, args.output_r1, args.output_r2, 
                          "TTTCTTATATGGGG", "TTGCTAGGACCGGCCTTAAAGC", log_file)

    log_message("Starting processing of matched reads...", log_file)
    process_reads(args.output_r1, args.output_r2, args.output_table, args.whitelist, log_file)

    log_message("Processing complete. Data saved to output files.", log_file)


if __name__ == "__main__":
    main()
