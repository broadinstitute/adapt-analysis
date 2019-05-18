"""Fetch accessions for taxonomic ID and a FASTA of sequences for those.
"""

import argparse

import ncbi
import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


INFLUENZA_TAXIDS = [11320, 11520]


def main(args):
    # Fetch accessions for tax_id
    segment = None if args.segment == 'None' else args.segment
    is_influenza = args.tax_id in INFLUENZA_TAXIDS
    accessions = ncbi.fetch_neighbors_acc(args.tax_id, segment,
            influenza=is_influenza)

    # Fetch fasta for these accessions (a tempfile)
    fasta_tf = ncbi.fetch_fastas(accessions)

    # Print out the contents of the FASTA
    with open(fasta_tf.name) as f:
        for line in f:
            print(line.rstrip())

    # Close the tempfile
    fasta_tf.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('tax_id', type=int,
            help="NCBI taxonomid ID")
    parser.add_argument('segment',
            help="Segment of taxonomy (or 'None' if unsegmented)")

    args = parser.parse_args()
    main(args)
