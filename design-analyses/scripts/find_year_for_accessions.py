"""Fetch accessions for taxonomic ID and group by year.
"""

import argparse
from collections import defaultdict
import re

import ncbi

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def fetch_collection_dates(accessions, batch_size=100, reqs_per_sec=2):
    """
    Args:
        accessions: collection of accessions

    Returns:
        dict {accession: [collection dates]}, where each value is
        a set in case there are multiple dates for an accession
    """
    xml_tf = ncbi.fetch_xml(accessions)
    source_features = ncbi.parse_genbank_xml_for_source_features(xml_tf.name)

    collection_dates = defaultdict(set)
    for accession, feats in source_features.items():
        for (name, value) in feats:
            if name == 'collection_date':
                collection_dates[accession].add(value)
    
    # Close the tempfile
    xml_tf.close()

    return collection_dates


def main(args):
    # Fetch accessions for tax_id
    segment = None if args.segment == 'None' else args.segment
    accessions = ncbi.fetch_neighbors_acc(args.tax_id, segment)

    # Fetch collection dates
    collection_dates = fetch_collection_dates(accessions)

    # Construct a pattern to match years in a date (1000--2999)
    year_p = re.compile('([1-2][0-9]{3})')

    # Parse a year from each collection date
    years = {}
    for acc in accessions:
        dates = collection_dates[acc]
        for date in dates:
            year_m = year_p.search(date)
            if year_m is None:
                # No year available; skip
                continue
            year = int(year_m.group(1))
            if acc in years and years[acc] != year:
                raise Exception(("Conflicting years for accession %s") % acc)
            years[acc] = year
        if acc not in years:
            # There is no date, or could not parse a year from the date
            years[acc] = 'Unknown'

    # Print a year for each accession (or 'Unknown')
    for acc in accessions:
        print('\t'.join([acc, str(years[acc])]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('tax_id', type=int)
    parser.add_argument('segment')

    args = parser.parse_args()
    main(args)
