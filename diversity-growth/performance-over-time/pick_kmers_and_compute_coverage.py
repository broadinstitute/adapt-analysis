"""Naively select a single k-mer (e.g., guide) from sequences
collected up to each year, and compute their coverage against
sequences from each year.
"""

import argparse
from collections import defaultdict
import statistics
import random
import re

from adapt import alignment
from adapt.prepare import ncbi_neighbors
from adapt.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def fetch_collection_dates(accessions, batch_size=100, reqs_per_sec=2):
    """Fetch sample collection dates for accessions.

    Args:
        accessions: collection of accessions

    Returns:
        dict {accession: [collection dates]}, where each value is
        a set in case there are multiple dates for an accession
    """
    # For accessions that are actually [accession].[version], strip
    # the .[version]
    accessions = [a.split('.')[0] for a in accessions]

    xml_tf = ncbi_neighbors.fetch_xml(accessions)
    source_features = ncbi_neighbors.parse_genbank_xml_for_source_features(xml_tf.name)

    collection_dates = defaultdict(set)
    for accession, feats in source_features.items():
        for (name, value) in feats:
            if name == 'collection_date':
                collection_dates[accession].add(value)
    
    # Close the tempfile
    xml_tf.close()

    return collection_dates


def parse_influenza_collection_dates(headers):
    """Parse collection date from influenza sequence header.

    The date is the third to last field, space-separated.

    Args:
        headers: collection of sequence headers

    Returns:
        dict {accession: str of collection date}
    """
    collection_dates = {}
    for header in headers:
        hs = header.split(' ')
        accver = hs[0]
        acc = accver.split('.')[0]
        date = hs[-3]
        collection_dates[acc] = date
    return collection_dates


def read_sequences(fn, influenza=True):
    """Read alignment and fetch collection dates.

    Args:
        fn: path to alignment of sequences; sequence headers
            contain accession
        influenza: if True, assume sequence headers are from NCBI influenza
            database

    Returns:
        tuple (Alignment object, dict {sequence index in alignment: year})
    """
    seqs = seq_io.read_fasta(fn)

    if influenza:
        collection_dates = parse_influenza_collection_dates(seqs.keys())
        # Below assumes a collection of possible dates, so transform the
        # dict
        collection_dates = {acc: [date] for acc, date in
                collection_dates.items()}
    else:
        accs = [h.split(' ')[0] for h in seqs.keys()]
        collection_dates = fetch_collection_dates(accs)

    seq_list = []
    i = 0
    years = {}
    seqs_skipped = 0
    for seq_header, seq in seqs.items():
        # Header may be [accession].[version] [details]
        acc = seq_header.split(' ')[0].split('.')[0]

        year = None
        for date in collection_dates[acc]:
            match = re.search('\d{4}', date)
            if match:
                year = int(match.group(0))
                break
        if year is None:
            # Cannot find year; skip this sequence
            seqs_skipped += 1
            continue

        seq_list += [seq]
        years[i] = year
        i += 1

    print(('Parsed year from %d sequences; skipped %d sequences') %
            (len(seq_list), seqs_skipped))

    aln = alignment.Alignment.from_list_of_seqs(seq_list)
    return (aln, years)


def seq_idxs_in_year_range(years, start, end=None):
    """Find indices of sequences whose year is within a range.

    Args:
        years: dict {sequence index: year}
        start/end: start (inclusive) and end (inclusove) year; if
            end is None, only pick ones from start year; if start
            is None, pick all sequences up to end

    Returns:
        set of indices from dates
    """
    idxs = set()
    for idx, year in years.items():
        if end is None:
            if year == start:
                idxs.add(idx)
        elif start is None:
            if year <= end:
                idxs.add(idx)
        else:
            if start <= year <= end:
                idxs.add(idx)
    return idxs


def bootstrap_alignment(aln):
    """Make new alignment representing bootstrap sample of sequences.

    This randomly samples with replacement aln.num_sequences from aln,
    and makes a new alignment object using those.

    Args:
        aln: alignment.Alignment object

    Returns:
        alignment.Alignment object
    """
    seq_list = aln.make_list_of_seqs()

    seq_list_resampled = random.choices(seq_list,
            k=aln.num_sequences)

    return alignment.Alignment.from_list_of_seqs(seq_list_resampled)


def find_mode_kmers(aln, k=30, num_kmers=15):
    """Find most common k-mers across alignment -- the ones found in
    the highest fraction of sequences.

    This requires that k-mers come from different sites (so only allows the
    most common one at each site) and, moreover, that they do not overlap.

    Args:
        aln: alignment.Alignment object
        k: k-mer length
        num_kmers: number of common k-mers to return

    Returns:
        list of k-mers that are most common across alignment
    """
    kmers_with_frac = []
    for pos in range(0, aln.seq_length - k + 1):
        aln_for_kmer = aln.extract_range(pos, pos + k)

        # When constructing k-mers, ignore any sequences in the alignment
        # that have a gap in this region
        seqs_with_gap = set(aln_for_kmer.seqs_with_gap())
        seqs_to_consider = set(range(aln_for_kmer.num_sequences)) - seqs_with_gap

        # Construct k-mer
        mode_kmer = aln_for_kmer.determine_most_common_sequence(
                seqs_to_consider=seqs_to_consider, skip_ambiguity=True)
        if mode_kmer is None:
            # No suitable k-mer; skip this position
            continue

        # Determine the fraction of the sequences that the k-mer is found in
        mode_kmer_bound = aln_for_kmer.sequences_bound_by_guide(
                mode_kmer, 0, 0, False)
        mode_kmer_frac = float(len(mode_kmer_bound)) / aln_for_kmer.num_sequences

        kmers_with_frac += [(mode_kmer, mode_kmer_frac, pos)]

    # Select the k-mers with the highest fraction bound
    selected_kmers = []
    best = sorted(kmers_with_frac, key=lambda x: x[1], reverse=True)
    for kmer, frac, pos in best:
        if len(selected_kmers) == num_kmers:
            break
        # Check if kmer overlaps with any in selected_kmers
        overlaps = False
        for other_kmer, other_kmer_pos in selected_kmers:
            if (pos <= other_kmer_pos + k and
                    other_kmer_pos <= pos + k):
                # kmer overlaps with other_kmer
                overlaps = True
                break
        if not overlaps:
            selected_kmers += [(kmer, pos)]
    return [x[0] for x in selected_kmers]


def find_mode_kmers_bootstrap(aln, k=30, sample_size=10):
    """Find mode k-mers over bootstrap samples of the alignment.

    Args:
        aln: alignment.Alignment object
        k: k-mer length
        sample_size: number of bootstrap samples

    Returns:
        list of sample_size lists of k-mers that are each mode across alignment
        for the bootstrap
    """
    mode_kmers_bootstrap = []
    for _ in range(sample_size):
        aln_bootstrap = bootstrap_alignment(aln)
        mode_kmers = find_mode_kmers(aln_bootstrap, k=k)
        mode_kmers_bootstrap += [mode_kmers]
    return mode_kmers_bootstrap


def extract_aln_with_seq_idxs(aln, seq_idxs):
    """Pull out alignment of sequences with given indices.

    Args:
        aln: alignment.Alignment object
        seq_idxs: set of indices of sequences in aln

    Returns:
        alignment.Alignment object
    """
    seq_list = aln.make_list_of_seqs(seqs_to_consider=seq_idxs)
    return alignment.Alignment.from_list_of_seqs(seq_list)


def make_designs_from_each_year(aln, years, start_year, end_year,
        included_years):
    """Produce "designs" for each year.

    This chooses k-mers (one per bootstrap sample) using all sequences
    collected recently for each year.

    Args:
        aln: alignment.Alignment object
        years: dict {sequence index: year}
        start_year/end_year: start year (inclusive) and end year
            (inclusive) to consider
        included_years: number of years to include in the design; for a
            year Y, this designs using sequence data from years
            (Y - included_years + 1) through Y

    Returns:
        dict {year: list of lists x_i of k-mers (one x_i per bootstrap)}
    """
    assert included_years >= 1

    # Choose k-mers for each year
    kmers_for_year = {}
    for year in range(start_year, end_year + 1):
        # Find all sequences collected in years up to and
        # including year
        print('Finding k-mers for year %d' % year)
        design_start_year = year - included_years + 1
        seq_idxs = seq_idxs_in_year_range(years, start=design_start_year, end=year)
        aln_year = extract_aln_with_seq_idxs(aln, seq_idxs)

        kmers_for_year[year] = find_mode_kmers_bootstrap(aln_year)
    return kmers_for_year


def compute_fraction_hit(aln, years, start_year, end_year, kmers, mismatches):
    """Find mean fraction of sequences hit by k-mers in each year.

    This tolerates mismatches between k-mer and target, but explicitly
    does not allow G-U pairing.

    The mean is calculated over the k-mers.

    Args:
        aln: alignment.Alignment object
        years: dict {sequence index: year}
        start_year/end_year: start year (inclusive) and end year
            (inclusive) to consider
        kmers: list of k-mers to query
        mismatches: number of mismatches to tolerate when determining
            if a k-mer 'hits' a target sequence

    Returns:
        dict {year: mean fraction hit of sequences in year hit}
    """
    k = len(kmers[0])

    frac_hit_in_year = {}
    for year in range(start_year, end_year + 1):
        # Get alignment containing only sequences from year
        seq_idxs = seq_idxs_in_year_range(years, start=year, end=None)
        aln_year = extract_aln_with_seq_idxs(aln, seq_idxs)

        # Scan along the alignment to find the site with
        # the most sequences bound for each k-mer (we don't have the site
        # where k-mer binds to, so we have to find it)
        fracs_bound = defaultdict(list)
        for pos in range(0, aln.seq_length - k + 1):
            aln_at_pos = aln_year.extract_range(pos, pos + k)

            for kmer in kmers:
                kmer_bound = aln_at_pos.sequences_bound_by_guide(
                        kmer, 0, mismatches, False)
                frac = float(len(kmer_bound)) / aln_year.num_sequences
                fracs_bound[kmer].append(frac)

        frac_bound_over_kmers = []
        for kmer in kmers:
            frac_bound = max(fracs_bound[kmer])
            frac_bound_over_kmers += [frac_bound]
        frac_hit_in_year[year] = statistics.mean(frac_bound_over_kmers)

    return frac_hit_in_year


def main(args):
    aln, years = read_sequences(args.in_fasta)
    start_year = min(years.values())
    end_year = max(years.values())

    # Produce designs
    if args.design_year:
        design_start_year, design_end_year = args.design_year, args.design_year
    else:
        design_start_year, design_end_year = start_year, end_year
    kmers_for_year = make_designs_from_each_year(aln, years,
            design_start_year, design_end_year, args.years_included_in_design)

    # Determine fraction of sequences from each year that
    # are hit by each k-mer, and write this to TSV file
    with open(args.out_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')

        header = ['design_year', 'sampling', 'kmers', 'test_year',
                'frac_hit']
        write_row(header)

        for design_year, kmers in kmers_for_year.items():
            for sample_i, sample_kmers in enumerate(kmers):
                frac_hit_in_year = compute_fraction_hit(
                        aln, years, start_year, end_year, sample_kmers,
                        args.mismatches)
                sample_kmers_str = ','.join(sample_kmers)
                for test_year, frac_hit in frac_hit_in_year.items():
                    row = [design_year, sample_i, sample_kmers_str,
                            test_year, frac_hit]
                    write_row(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--in-fasta',
        required=True,
        help=("Alignment of sequences"))
    parser.add_argument('--out-tsv',
        required=True,
        help=("Path to output TSV"))
    parser.add_argument('--mismatches',
        type=int,
        default=1,
        help=("Number of mismatches when determining whether a k-mer "
              "hits/detects a target sequence"))
    parser.add_argument('--design-year',
        type=int,
        help=("If set, only produce and evaluate designs produced in this "
              "year"))
    parser.add_argument('--years-included-in-design',
        type=int,
        default=2,
        help=("Number of years to include in a design. When designing for "
              "a year Y, this uses sequence data from all years in [Y - "
              "YEARS_INCLUDED_IN_DESIGN + 1, Y]."))

    args = parser.parse_args()
    main(args)
