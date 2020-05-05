"""Compute entropy of k-mers in sliding window at different time points.

This uses the alignment module from ADAPT, so the ADAPT environment must be loaded.

This is written specific to parsing the format of GISAID sequence headers for
SARS-CoV-2.
"""

import argparse
import datetime
import random

import numpy as np

from adapt import alignment
from adapt.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def read_sequences(fn, ref_acc):
    """Read alignment and dates.

    Args:
        fn: path to alignment of sequences from GISAID; sequence headers contain
            sample date
        ref_acc: accession of reference sequence found in header

    Returns:
        tuple (Alignment object, dict {sequence index in alignment: date},
        index of reference sequence)
    """
    seqs = seq_io.read_fasta(fn)

    seq_list = []
    i = 0
    dates = {}
    seqs_skipped = 0
    ref_idx = None
    for seq_header, seq in seqs.items():
        # Last field, split by '|', should have date
        date_str = seq_header.split('|')[-1]
        try:
            date = datetime.datetime.strptime(date_str, '%Y-%m-%d')
        except ValueError:
            # The string may only contain month (not day); we could
            # try parsing that. For now, just skip the sequence
            seqs_skipped += 1
            continue

        if ref_acc in seq_header:
            assert ref_idx is None
            ref_idx = i

        seq_list += [seq]
        dates[i] = date
        i += 1

    print(('Parsed date from %d sequences; skipped %d sequences') %
            (len(seq_list), seqs_skipped))
    assert ref_idx is not None

    aln = alignment.Alignment.from_list_of_seqs(seq_list)
    return (aln, dates, ref_idx)


def ref_aln_pos_map(aln, ref_idx):
    """Make map of position in reference sequence to alignment position.

    Args:
        aln: alignment.Alignment object
        ref_idx: index of reference sequence in aln

    Returns:
        list x such that x[i] gives the position in the alignment of the
        i'th position in the reference sequence
    """
    # Get the reference sequence
    s = aln.make_list_of_seqs(seqs_to_consider={ref_idx},
            remove_gaps=False)[0]
    s_without_gaps = s.replace('-', '')

    m = [-1 for _ in range(len(s_without_gaps))]
    curr_s_pos = 0
    for i in range(len(s)):
        if s[i] == '-':
            # gap
            continue
        else:
            # nucleotide in s
            m[curr_s_pos] = i
            curr_s_pos += 1
    assert curr_s_pos == len(s_without_gaps)
    return m


def seqs_idx_in_date_range(dates, start, end):
    """Find indices of sequences whose date is within a range.

    Args:
        dates: dict {sequence index: date}
        start/end: start (inclusive) and end (exclusive) dates

    Returns:
        set of indices from dates
    """
    idxs = set()
    for idx, date in dates.items():
        if start <= date < end:
            idxs.add(idx)
    return idxs


def entropy_of_kmers(kmers, base=2.0):
    """Compute entropy of a list of k-mers.

    Args:
        kmers: list of k-mers
        base: base of logarithm

    Returns:
        float giving entropy
    """
    unique, counts = np.unique(kmers, return_counts=True)
    frequencies = counts / counts.sum()
    return np.sum(-1.0 * frequencies * np.log(frequencies)/np.log(base))


def mean_entropy_in_window(ref_start, ref_end, ref_pos_map,
        aln, seq_idxs, k=28):
    """Find mean entropy across sites in a window.

    Args:
        ref_start/ref_end: start (inclusive) and end (exclusive)
            positions of a window in the reference sequence
        ref_pos_map: dict {position in reference: position in aln}
        aln: alignment.Alignment object
        seq_idxs: indices of sequences to use
        k: k-mer length

    Returns:
        float giving entropy
    """
    entropies = []
    for pos in range(ref_start, ref_end - k):
        # Find the corresponding range in the alignment
        aln_start, aln_end = ref_pos_map[pos], ref_pos_map[pos + k]
        # Extract site of k-mer
        extract = aln.extract_range(aln_start, aln_end)

        # Get all sequences at this site
        kmers = extract.make_list_of_seqs(seqs_to_consider=seq_idxs,
                remove_gaps=True)

        # Compute entropy
        entropy = entropy_of_kmers(kmers)
        entropies += [entropy]
    return np.mean(entropies)


def main(args):
    aln, dates, ref_idx = read_sequences(args.in_fasta, args.ref_acc)
    ref_pos_map = ref_aln_pos_map(aln, ref_idx)
    ref_len = len(ref_pos_map)

    start_date = datetime.datetime.strptime(args.start_date, '%Y-%m-%d')
    end_date = datetime.datetime.strptime(args.end_date, '%Y-%m-%d')
    date_incr = datetime.timedelta(days=args.interval_days)

    rows = []

    # Iterate over date ranges
    date = start_date
    while date + date_incr <= end_date:
        date_range_start = date
        date_range_end = date + date_incr

        date_range_start_str = date_range_start.strftime('%Y-%m-%d')
        date_range_end_str = date_range_end.strftime('%Y-%m-%d')

        print("On date range: [%s, %s)" % (date_range_start_str,
            date_range_end_str))

        # Find sequences in this date range
        seq_idxs = seqs_idx_in_date_range(dates, date_range_start,
                date_range_end)
        if args.sample_size_per_date_interval > len(seq_idxs):
            # Too few sequences to sample; warn and move on
            print(("Number of sequences in date range [%s, %s) is "
                "%d, which is fewer than the sample size (%d); not sampling") %
                (date_range_start_str, date_range_end_str, len(seq_idxs),
                    args.sample_size_per_date_interval))
        else:
            # Sample randomly (without replacement)
            seq_idxs = random.choices(list(seq_idxs),
                    k=args.sample_size_per_date_interval)
            seq_idxs = set(seq_idxs)

        # Iterate over genome windows
        i = 0
        while i + args.window_length < ref_len:
            window_start = i
            window_end = i + args.window_length

            # Get the mean entropy in this window
            entropy = mean_entropy_in_window(window_start, window_end,
                    ref_pos_map, aln, seq_idxs)

            # Record the entropy
            rows += [(date_range_start_str, date_range_end_str, window_start,
                window_end, entropy)]

            # Increment the window start by the stride
            i += args.window_stride

        # Increment start of date range
        date += date_incr

    # Write the output
    with open(args.out_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')

        header = ['start_date', 'end_date', 'window_start', 'window_end',
                'entropy']
        write_row(header)

        for row in rows:
            write_row(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--in-fasta',
        required=True,
        help=("Path to FASTA of aligned sequences from GISAID"))
    parser.add_argument('--out-tsv',
        required=True,
        help=("Path to output TSV to which to write results"))
    parser.add_argument('--ref-acc',
        required=True,
        help=("Accession of reference; positions are relative to this"))
    parser.add_argument('--start-date',
        required=True,
        help=("Start date in YYYY-MM-DD format"))
    parser.add_argument('--end-date',
        required=True,
        help=("End date in YYYY-MM-DD format"))
    parser.add_argument('--interval-days',
        default=7,
        type=int,
        help=("Number of days in a time window"))
    parser.add_argument('--window-length',
        default=500,
        type=int,
        help=("Window length, in nucleotides, to average over"))
    parser.add_argument('--window-stride',
        default=500,
        type=int,
        help=("Amoutn by which to stride window, in nucleotides"))
    parser.add_argument('--sample-size-per-date-interval',
        default=100,
        type=int,
        help=("Number of sequences to sample randomly from each "
              "date interval, to remove temporal sampling biases"))

    args = parser.parse_args()
    main(args)
