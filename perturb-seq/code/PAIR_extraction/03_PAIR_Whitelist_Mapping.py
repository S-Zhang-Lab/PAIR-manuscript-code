import argparse
import csv
from Levenshtein import distance
from collections import defaultdict
import time

def create_kmer_index(whitelist_pairs, k):
    """Create a k-mer index for faster matching."""
    kmer_index = defaultdict(list)
    for pair_seq, data in whitelist_pairs.items():
        # Generate k-mers for each whitelist sequence
        for i in range(len(pair_seq) - k + 1):
            kmer = pair_seq[i:i + k]
            if 'N' not in kmer:  # Only index k-mers without N
                kmer_index[kmer].append((pair_seq, data))
    return kmer_index

def efficient_kmer_match(input_string, kmer_index, whitelist_pairs, k, tolerance):
    """Improved k-mer matching using index lookup."""
    best_match = None
    best_distance = float('inf')
    
    # Get candidate matches using k-mer index
    candidates = set()
    for i in range(len(input_string) - k + 1):
        kmer = input_string[i:i + k]
        if 'N' not in kmer and kmer in kmer_index:
            for pair_seq, data in kmer_index[kmer]:
                candidates.add(pair_seq)
    
    # Check candidates with Levenshtein distance
    for candidate in candidates:
        edit_dist = distance(input_string, candidate)
        if edit_dist <= tolerance and edit_dist < best_distance:
            best_distance = edit_dist
            best_match = whitelist_pairs[candidate]
            
            # Early exit if exact match found
            if edit_dist == 0:
                break
    
    return best_match, best_distance

def match_pairs(input_file, whitelist_file, output_file, k, tolerance):
    """Improved mapping function with k-mer indexing."""
    # Load whitelist pairs with concatenated sequences
    print(f"Loading whitelist file: {whitelist_file}")
    whitelist_pairs = {}
    with open(whitelist_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            # Create concatenated sequence from CRISPRa and CasRx
            pair_seq = row['CRISPRa'] + row['CasRx']
            whitelist_pairs[pair_seq] = row
    print(f"Whitelist loaded. Total entries: {len(whitelist_pairs)}")
    
    # Create k-mer index
    print(f"Creating k-mer index with k={k}...")
    start_time = time.time()
    kmer_index = create_kmer_index(whitelist_pairs, k)
    print(f"K-mer index created in {time.time() - start_time:.2f} seconds")
    
    # Process input file
    total_rows = 0
    mapped_rows = 0
    start_time = time.time()
    
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        fieldnames = reader.fieldnames + ['Mapped_PAIR_seq', 'Mapped_name', 'Edit_Distance']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        for row in reader:
            total_rows += 1
            # Concatenate input sequences
            pair = row['CRISPRa_gRNA'] + row['CasRx_crRNA']
            match_data, min_distance = efficient_kmer_match(pair, kmer_index, whitelist_pairs, k, tolerance)
            
            if match_data:
                row['Mapped_PAIR_seq'] = match_data['CRISPRa'] + match_data['CasRx']
                row['Mapped_name'] = f"{match_data['CRISPRa_name']}_{match_data['CasRx_name']}"
                row['Edit_Distance'] = min_distance
                writer.writerow(row)
                mapped_rows += 1
            
            if total_rows % 1000 == 0:
                elapsed_time = time.time() - start_time
                rate = total_rows / elapsed_time
                print(f"Processed {total_rows} rows ({rate:.1f} rows/sec). Mapped: {mapped_rows}")
    
    print(f"Finished processing. Total: {total_rows}, Mapped: {mapped_rows}")
    print(f"Total time: {time.time() - start_time:.2f} seconds")
    return total_rows, mapped_rows

def main():
    parser = argparse.ArgumentParser(description="Optimized PAIR mapping using k-mer indexing")
    parser.add_argument('--input', required=True, help="Input data file")
    parser.add_argument('--whitelist', required=True, help="Whitelist file")
    parser.add_argument('--output', required=True, help="Output file")
    parser.add_argument('--k', type=int, default=10, help="K-mer size")
    parser.add_argument('--tolerance', type=int, default=1, help="Maximum allowed Levenshtein distance")
    parser.add_argument('--log', required=False, help="Log file")
    args = parser.parse_args()

    # Direct mapping without intermediate files
    total_rows, mapped_rows = match_pairs(args.input, args.whitelist, args.output, args.k, args.tolerance)

    # Log statistics
    stats = {
        "Total rows": total_rows,
        "Mapped rows": mapped_rows,
        "Mapping rate": f"{(mapped_rows/total_rows*100):.2f}%" if total_rows > 0 else "0%",
        "K-mer size": args.k,
        "Tolerance": args.tolerance
    }
    
    if args.log:
        with open(args.log, 'w') as log_file:
            for key, value in stats.items():
                log_file.write(f"{key}: {value}\n")

if __name__ == "__main__":
    main()