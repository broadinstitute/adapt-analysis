"""Fetch a FASTA of sequences for accessions.
"""

import argparse

from adapt.prepare import ncbi_neighbors

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    # Read accessions
    accessions = []
    with open(args.accessions) as f:
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                continue
            accessions += [line]

    # Fetch fasta for these accessions (a tempfile)
    fasta_tf = ncbi_neighbors.fetch_fastas(accessions)

    # Print out the contents of the FASTA
    with open(fasta_tf.name) as f:
        for line in f:
            print(line.rstrip())

    # Close the tempfile
    fasta_tf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('accessions',
            help="Path to file listing accessions")

    args = parser.parse_args()
    main(args)
