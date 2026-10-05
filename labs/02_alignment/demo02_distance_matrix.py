"""Nucleotide p-distance from global pairwise alignments or a supplied MSA.

Raw input: --fasta data/sample/tp53_dna_multi.fasta
Aligned FASTA: --fasta your_msa.fasta --aligned
Exclude columns with gaps/ambiguous bases (pairwise deletion).
This is an observed mismatch fraction, not a corrected evolutionary distance.
"""
import argparse
import csv
from itertools import combinations
import sys

from Bio import Align, SeqIO


def aligned_counts(a, b):
    if len(a) != len(b):
        raise ValueError("Aligned sequences must have equal lengths; truncation is not alignment.")
    pairs = [(x, y) for x, y in zip(a.upper(), b.upper()) if x in "ACGT" and y in "ACGT"]
    return sum(x != y for x, y in pairs), len(pairs)


def pairwise_counts(a, b):
    aligner = Align.PairwiseAligner(mode="global", match_score=1, mismatch_score=-1,
                                   open_gap_score=-2, extend_gap_score=-2)
    alignment = aligner.align(a, b)[0]
    mismatches = compared = 0
    for (a0, a1), (b0, b1) in zip(*alignment.aligned):
        diff, count = aligned_counts(a[a0:a1], b[b0:b1])
        mismatches += diff
        compared += count
    return mismatches, compared


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--aligned", action="store_true", help="Input is an existing aligned FASTA")
    args = parser.parse_args()
    with open(args.fasta) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    if len(records) < 2 or any(not len(record) for record in records):
        parser.error("Provide at least two nonempty nucleotide sequences.")
    if len({r.id for r in records}) != len(records):
        parser.error("Sequence IDs must be unique.")
    seqs = [str(r.seq).upper() for r in records]
    allowed = set("ACGTRYSWKMBDHVN" + ("-" if args.aligned else ""))
    if any(set(seq) - allowed for seq in seqs):
        parser.error("Expected nucleotide IUPAC symbols (T convention); gaps require --aligned.")
    if args.aligned and len({len(s) for s in seqs}) != 1:
        parser.error("All MSA rows must have equal lengths.")
    if not args.aligned and max(map(len, seqs)) > 10000:
        parser.error("Demo limit: 10,000 nt per raw sequence. Use shorter comparable regions or an MSA.")
    mode = "provided MSA" if args.aligned else "global pairwise; match=1 mismatch=-1 gap=-2"
    print(f"# {mode}; exclude gaps/ambiguous bases; denominator varies by pair", file=sys.stderr)
    writer = csv.writer(sys.stdout)
    writer.writerow(["id1", "id2", "mismatches", "p_distance", "compared_sites"])
    for i, j in combinations(range(len(seqs)), 2):
        method = aligned_counts if args.aligned else pairwise_counts
        differences, count = method(seqs[i], seqs[j])
        distance = f"{differences / count:.4f}" if count else "NA"
        writer.writerow([records[i].id, records[j].id, differences, distance, count])


if __name__ == "__main__":
    main()
