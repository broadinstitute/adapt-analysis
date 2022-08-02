"""Make table of sequences over time, given an accession list.

This uses entry creation dates for the year, rather than
sample collection dates, in case there is bias in sample
collection dates (e.g., people did not start using them
until more recent); not every entry has a collection date,
but every entry on GenBank has a entry creation date.

Note that this counts any sequence as its own -- i.e., any
chromosome or segment gets counted as one.
"""

import argparse
from collections import Counter
from collections import defaultdict
from os import environ
import requests

from adapt.prepare import align
from adapt.prepare import ncbi_neighbors

__author__ = 'Hayden Metsky <hayden@mit.edu>'


INCLUDE_KMER_COUNTS = False

START_YEAR = 1990
END_YEAR = 2022


# Set NCBI api key to environment variable NCBI_API_KEY
ncbi_api_key = environ.get('NCBI_API_KEY')
if ncbi_api_key is not None:
    print("Setting NCBI API key from environment variable")
    ncbi_neighbors.set_ncbi_api_key(ncbi_api_key)
else:
    print("NCBI API key is not set in environment variable")


def num_unique_kmers(accessions, k=31):
    """Download sequences and count number of unique k-mers.

    This uses k=31 by default, which is what kraken uses for metagenomic
    classification.

    Args:
        accessions: list of accessions
        k: k-mer size

    Returns:
        number of unique k-mers
    """
    accessions = list(accessions)
    seqs_fp = ncbi_neighbors.fetch_fastas(accessions)
    seqs = align.read_unaligned_seqs(seqs_fp)
    all_kmers = set()
    for accver, seq in seqs.items():
        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        all_kmers |= set(kmers)
    seqs_fp.close()

    return len(all_kmers)


def main(args):
    # Read accessions
    print("Reading accessions")
    accessions = []
    with open(args.accession_list) as f:
        for line in f:
            # Strip out a version number (after a period), if one is present
            acc = line.rstrip().split('.')[0]
            accessions += [acc]

    # Fetch metadata
    print("Fetching metadata")
    metadata = ncbi_neighbors.fetch_metadata(accessions)

    # Fetch sequences and compute number of unique k-mers
    print("Determining number of unique k-mers")
    if INCLUDE_KMER_COUNTS:
        # Do this for the collection of accessions from each year
        acc_for_year = defaultdict(set)
        for acc, m in metadata.items():
            acc_year = m['entry_create_year']
            acc_for_year[acc_year].add(acc)
        nuk_for_year = defaultdict(int)
        for year, accessions in acc_for_year.items():
            nuk = num_unique_kmers(accessions)
            nuk_for_year[year] = nuk

    # Write TSV
    print("Writing output TSVs")
    with open(args.out_per_year_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')

        header = ['year', 'num_sequences']
        if INCLUDE_KMER_COUNTS:
            header += ['num_unique_kmers']
        write_row(header)

        year_sequence_count = defaultdict(int)
        for acc, m in metadata.items():
            year_sequence_count[m['entry_create_year']] += 1
        for year in range(START_YEAR, END_YEAR+1):
            row = [year, year_sequence_count[year]]
            if INCLUDE_KMER_COUNTS:
                row += [nuk_for_year[year]]
            write_row(row)

    with open(args.out_cumulative_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')

        header = ['year', 'cumulative_num_sequences']
        if INCLUDE_KMER_COUNTS:
            header += ['cumulative_num_unique_kmers']
        write_row(header)

        for year in range(START_YEAR, END_YEAR+1):
            year_sequence_count = 0
            for acc, m in metadata.items():
                if m['entry_create_year'] <= year:
                    year_sequence_count += 1
            row = [year, year_sequence_count]
            if INCLUDE_KMER_COUNTS:
                year_kmer_count = 0
                for y, nuk in num_for_year.items():
                    if y <= year:
                        year_kmer_count += nuk
                row += [year_kmer_count]
            write_row(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('accession_list')
    parser.add_argument('out_per_year_tsv')
    parser.add_argument('out_cumulative_tsv')
    args = parser.parse_args()
    main(args)
