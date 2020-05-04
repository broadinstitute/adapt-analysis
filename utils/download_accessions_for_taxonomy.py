"""Fetch accessions for taxonomic ID and a FASTA of sequences for those.
"""

import argparse

import ncbi
import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


INFLUENZA_TAXIDS = [11320, 11520, 11552]


def main(args):
    # Fetch accessions for tax_id
    segment = None if args.segment == 'None' else args.segment
    is_influenza = args.tax_id in INFLUENZA_TAXIDS
    accessions = ncbi.fetch_neighbors_acc(args.tax_id, segment,
            influenza=is_influenza)

    if len(accessions) == 0:
        raise Exception("No accessions for (%d, %s)" % (args.tax_id,
            args.segment))

    for acc in accessions:
        row = [args.tax_id, args.segment, acc]
        print('\t'.join(str(x) for x in row))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('tax_id', type=int,
            help="NCBI taxonomid ID")
    parser.add_argument('segment',
            help="Segment of taxonomy (or 'None' if unsegmented)")

    args = parser.parse_args()
    main(args)
