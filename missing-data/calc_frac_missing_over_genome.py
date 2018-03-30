"""At each site in an alignment, calculate the fraction of sequences with
missing data.

This assumes 'N' to be missing data, as well as '-' but only on the ends.
It assumes these to be the only symbols of missing data.
"""

import argparse
import sys

sys.path.append("../utils")
import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    aln = seq_io.read_fasta(args.in_fasta)

    # Verify this is an alignment (all sequences are the same
    # length) and get the length
    seq_len = None
    for name, seq in aln.items():
        if seq_len is None:
            seq_len = len(seq)
        assert seq_len == len(seq)

    # Count the number of sequences at each position that have
    # missing data
    missing_count = [0 for _ in range(seq_len)]
    for name, seq in aln.items():
        # Find the start position (everything before this will
        # be counted as missing data)
        start_pos = 0
        for i in range(seq_len):
            if seq[i] != '-':
                start_pos = i
                break

        # Find the end position (everything after this will be
        # counted as missing data)
        end_pos = seq_len - 1
        for i in reversed(range(seq_len)):
            if seq[i] != '-':
                end_pos = i
                break

        for i in range(seq_len):
            if i < start_pos or i > end_pos or seq[i] == 'N':
                missing_count[i] += 1


    num_seqs = len(aln)
    missing_frac = [float(missing_count[i]) / num_seqs for i in range(seq_len)]

    # Print the fraction missing at each position
    header = ['pos', 'frac.missing']
    with open(args.out_tsv, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for i in range(seq_len):
            row = [str(x) for x in [i+1, missing_frac[i]]]
            f.write('\t'.join(row) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('in_fasta', help=("Alignment in fasta format"))
    parser.add_argument('out_tsv', help=("Output table in TSV format"))

    args = parser.parse_args()
    main(args)
