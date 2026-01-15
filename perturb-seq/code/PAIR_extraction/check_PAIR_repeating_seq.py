#!/usr/bin/env python3
"""
Scan the first N reads of a gzipped FASTQ and look for occurrences of two
whitelist sequences allowing up to K mismatches (default K=1).

Outputs:
  1) A per-read table of match positions (0-based and 1-based) for each pattern.
  2) A summary JSON with frequencies, position histograms, and mismatch histograms.

Usage:
  python scan_whitelist_fastq_mm.py \
    --fastq PAIR_TT_C12-SI_TT_C12_2333FLLT3_S3_L007_R2_001.fastq.gz \
    --n_reads 100 \
    --max_mismatch 1 \
    --out_prefix whitelist_scan
"""

import argparse
import gzip
import json
from collections import Counter
from dataclasses import dataclass
from typing import Dict, List, Tuple, Iterator


def iter_fastq_gz(path: str) -> Iterator[Tuple[str, str, str, str]]:
    """Yield FASTQ records (header, seq, plus, qual) from a .fastq.gz file."""
    with gzip.open(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                return
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            if not qual:
                return
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), qual.rstrip("\n")


def mismatches_leq_k(a: str, b: str, k: int) -> int:
    """
    Return mismatch count if <=k else (k+1). Early exit for speed.
    Assumes len(a) == len(b).
    """
    mm = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            mm += 1
            if mm > k:
                return k + 1
    return mm


def find_all_positions_upto_k_mismatch(seq: str, pat: str, k: int) -> Tuple[List[int], List[int]]:
    """
    Find all (possibly overlapping) occurrences of pat in seq allowing <=k mismatches.
    Returns:
      - positions_0based: list of start indices
      - mismatches: list of mismatch counts aligned with positions
    Optimized with a seed-and-verify approach:
      - use exact matches of a long seed (half the pattern) to generate candidates
      - verify full pattern with early-exit mismatch counter
    """
    n = len(seq)
    m = len(pat)
    if m == 0 or n < m:
        return [], []

    # Choose a seed: take the longer half (works well for k=1)
    # We use both halves to reduce false negatives when mismatch is inside the seed.
    half = m // 2
    seed1 = pat[: m - half]      # longer half if m is odd
    seed2 = pat[m - half :]      # shorter half

    # offsets of the seeds within the pattern
    off1 = 0
    off2 = m - half

    candidates = set()

    # Find all occurrences of seed1, translate to candidate start positions
    start = 0
    while True:
        idx = seq.find(seed1, start)
        if idx == -1:
            break
        cand = idx - off1
        if 0 <= cand <= n - m:
            candidates.add(cand)
        start = idx + 1

    # Find all occurrences of seed2, translate to candidate start positions
    start = 0
    while True:
        idx = seq.find(seed2, start)
        if idx == -1:
            break
        cand = idx - off2
        if 0 <= cand <= n - m:
            candidates.add(cand)
        start = idx + 1

    # If no candidates from seeds, we can bail early
    if not candidates:
        return [], []

    # Verify candidates
    positions = []
    mism_counts = []
    for cand in sorted(candidates):
        window = seq[cand : cand + m]
        mm = mismatches_leq_k(window, pat, k)
        if mm <= k:
            positions.append(cand)
            mism_counts.append(mm)

    return positions, mism_counts


@dataclass
class PatternStats:
    pattern: str
    reads_with_match: int = 0
    total_matches: int = 0

    # histograms
    start_pos_0based: Counter = None
    start_pos_1based: Counter = None
    matches_per_read: Counter = None
    mismatches_hist: Counter = None  # how many matches had 0 vs 1 mismatch, etc.

    def __post_init__(self):
        self.start_pos_0based = Counter()
        self.start_pos_1based = Counter()
        self.matches_per_read = Counter()
        self.mismatches_hist = Counter()

    def update(self, positions_0based: List[int], mismatches: List[int]) -> None:
        assert len(positions_0based) == len(mismatches)

        if positions_0based:
            self.reads_with_match += 1
            self.matches_per_read[len(positions_0based)] += 1
            self.total_matches += len(positions_0based)

            for p0, mm in zip(positions_0based, mismatches):
                self.start_pos_0based[p0] += 1
                self.start_pos_1based[p0 + 1] += 1
                self.mismatches_hist[mm] += 1
        else:
            self.matches_per_read[0] += 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True, help="Path to .fastq.gz file")
    ap.add_argument(
        "--n_reads",
        type=int,
        default=100,
        help="Number of reads to scan (0 = all reads; default: 100)",
    )
    ap.add_argument("--out_prefix", default="whitelist_scan", help="Output prefix")
    ap.add_argument("--max_mismatch", type=int, default=1, help="Max mismatches allowed (default: 1)")
    ap.add_argument(
        "--seq1",
        default="TTGCTAGGACCGGCCTTAAAGC",
        help="Whitelist sequence 1",
    )
    ap.add_argument(
        "--seq2",
        default="GTTTCAAACCCCGACCAGTTGGTAGGGGTTTACTTG",
        help="Whitelist sequence 2",
    )
    args = ap.parse_args()

    k = args.max_mismatch
    if k < 0:
        raise ValueError("--max_mismatch must be >= 0")

    patterns = {"SEQ1": args.seq1.upper(), "SEQ2": args.seq2.upper()}
    stats: Dict[str, PatternStats] = {lab: PatternStats(pat) for lab, pat in patterns.items()}

    per_read_rows = []
    scanned = 0

    for header, seq, plus, qual in iter_fastq_gz(args.fastq):
        if args.n_reads and scanned >= args.n_reads:
            break

        seq_u = seq.upper()
        read_id = header.split()[0].lstrip("@")

        row = {
            "read_index": scanned + 1,  # 1-based read index
            "read_id": read_id,
            "seq_len": len(seq_u),
        }

        for label, pat in patterns.items():
            pos0, mm = find_all_positions_upto_k_mismatch(seq_u, pat, k)
            stats[label].update(pos0, mm)

            row[f"{label}_n_matches"] = len(pos0)
            row[f"{label}_pos0"] = ",".join(map(str, pos0)) if pos0 else ""
            row[f"{label}_pos1"] = ",".join(str(p + 1) for p in pos0) if pos0 else ""
            row[f"{label}_mm"] = ",".join(map(str, mm)) if mm else ""  # mismatch counts aligned with positions

        per_read_rows.append(row)
        scanned += 1

    # Write per-read TSV
    tsv_path = f"{args.out_prefix}.per_read.tsv"
    cols = [
        "read_index", "read_id", "seq_len",
        "SEQ1_n_matches", "SEQ1_pos0", "SEQ1_pos1", "SEQ1_mm",
        "SEQ2_n_matches", "SEQ2_pos0", "SEQ2_pos1", "SEQ2_mm",
    ]
    with open(tsv_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in per_read_rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # Summary JSON
    summary = {
        "fastq": args.fastq,
        "n_reads_requested": args.n_reads,
        "n_reads_scanned": scanned,
        "max_mismatch": k,
        "patterns": {},
    }

    for label, st in stats.items():
        summary["patterns"][label] = {
            "pattern": st.pattern,
            "reads_with_match": st.reads_with_match,
            "reads_scanned": scanned,
            "freq_reads_with_match": (st.reads_with_match / scanned) if scanned else 0.0,
            "total_matches": st.total_matches,
            "avg_matches_per_read": (st.total_matches / scanned) if scanned else 0.0,
            "start_pos_hist_0based": dict(sorted(st.start_pos_0based.items())),
            "start_pos_hist_1based": dict(sorted(st.start_pos_1based.items())),
            "matches_per_read_hist": dict(sorted(st.matches_per_read.items())),
            "mismatches_hist": dict(sorted(st.mismatches_hist.items())),  # e.g., {0:..., 1:...}
        }

    json_path = f"{args.out_prefix}.summary.json"
    with open(json_path, "w") as out:
        json.dump(summary, out, indent=2)

    print(f"Done. Scanned {scanned} reads.")
    print(f"Wrote per-read positions: {tsv_path}")
    print(f"Wrote summary stats:     {json_path}")


if __name__ == "__main__":
    main()