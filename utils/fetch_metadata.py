"""Download metadata for a list of accessions.
"""

import argparse

import ncbi


def make_table(metadata, out_tsv):
    with open(out_tsv, 'w') as fw:
        header = ['accession', 'country', 'year', 'entry_create_year',
                'taxid']
        fw.write('\t'.join(header) + '\n')
        for accession in metadata:
            row = [accession] + [str(metadata[accession][c]) for c in header[1:]]
            fw.write('\t'.join(row) + '\n')


def main(args):
    metadata = ncbi.fetch_metadata(args.a)
    make_table(metadata, args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', nargs='+',
            help=("Space-separated list of accessions"))
    parser.add_argument('-o', help=("Output TSV"))
    args = parser.parse_args()

    main(args)
