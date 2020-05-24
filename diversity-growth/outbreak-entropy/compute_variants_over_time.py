"""Compute variants across genome at different time points.

This uses the alignment module from ADAPT, so the ADAPT environment must be loaded.

This is written specific to parsing the format of GISAID sequence headers for
SARS-CoV-2.
"""

import argparse
import datetime
import random

import numpy as np

from compute_entropy_over_time import read_sequences
from compute_entropy_over_time import ref_aln_pos_map
from compute_entropy_over_time import seqs_idx_in_date_range
from adapt import alignment
from adapt.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def find_variants(aln, ref_idx, ref_pos_map):
    """Find variants against a reference genome.

    Args:
        aln: alignment.Alignment object
        ref_idx: index of reference in aln
        ref_pos_map: list x where x[i] is position in aln of i'th position
            in the reference

    Returns:
        list x where x[i] is a dict {sequence index in aln: False/True
        indicating whether sequence index has a variant (True) against
        the reference at position i of the reference}
    """
    ref_seq = aln.make_list_of_seqs(seqs_to_consider={ref_idx})[0]

    variants = []
    for ref_pos, aln_pos in enumerate(ref_pos_map):
        ref_base = ref_seq[ref_pos]
        has_variant = {}

        # Get a string of all the bases at aln_pos, and iterate over
        # this
        aln_at_pos = aln.seqs[aln_pos]
        for seq_idx, b in enumerate(aln_at_pos):
            if b.upper() not in ['A', 'C', 'G', 'T']:
                # b may be a gap or ambiguity; skip it
                continue
            has_variant[seq_idx] = (b != ref_base)

        variants += [has_variant]

    return variants


def main(args):
    aln, dates, ref_idx = read_sequences(args.in_fasta, args.ref_acc)
    ref_pos_map = ref_aln_pos_map(aln, ref_idx)
    ref_len = len(ref_pos_map)

    variants = find_variants(aln, ref_idx, ref_pos_map)

    start_date = datetime.datetime.strptime(args.start_date, '%Y-%m-%d')
    end_date = datetime.datetime.strptime(args.end_date, '%Y-%m-%d')
    date_incr = datetime.timedelta(days=args.interval_days)

    rows = []

    # Iterate over date ranges
    date = start_date
    while date + date_incr <= end_date:
        date_str = date.strftime('%Y-%m-%d')

        print("On date: %s" % date_str)

        if args.cumulative:
            # Use everything up to this date
            date_range_start = datetime.datetime.min
            date_range_end = date
        else:
            # Use everything between the previous date and this one
            if date == start_date:
                # There is no previous date; skip this interval
                date += date_incr
                continue
            date_range_start = date - date_incr
            date_range_end = date

        # Find sequences in this date range
        seq_idxs = seqs_idx_in_date_range(dates, date_range_start,
                date_range_end)
        if args.sample_size_per_date_interval <= 0:
            # Do not subsample; use all
            seq_idxs = set(seq_idxs)
        elif args.sample_size_per_date_interval > len(seq_idxs):
            # Too few sequences to sample; warn and move on
            date_range_start_str = date_range_start.strftime('%Y-%m-%d')
            date_range_end_str = date_range_end.strftime('%Y-%m-%d')
            print(("Number of sequences in date range [%s, %s) is "
                "%d, which is fewer than the sample size (%d); not sampling") %
                (date_range_start_str, date_range_end_str, len(seq_idxs),
                    args.sample_size_per_date_interval))
        else:
            # Sample randomly (without replacement)
            seq_idxs = random.sample(list(seq_idxs),
                    args.sample_size_per_date_interval)
            seq_idxs = set(seq_idxs)

        # Count the number of variants against the reference at every
        # position of the reference
        for pos, d in enumerate(variants):
            num_counted = 0
            num_with_variants = 0
            for seq_idx, has_variant in d.items():
                if seq_idx not in seq_idxs:
                    # seq_idx is not sampled; ignore it
                    continue
                num_counted += 1
                if has_variant:
                    num_with_variants += 1
            rows += [(date_str, pos, num_with_variants, num_counted)]

        # Increment start of date range
        date += date_incr

    # Write the output
    with open(args.out_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')

        header = ['date', 'pos', 'num_genomes_with_variant',
                  'num_genomes_counted']
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
        help=("Accession of reference; positions and variants are relative to "
              "this"))
    parser.add_argument('--start-date',
        required=True,
        help=("Start date in YYYY-MM-DD format"))
    parser.add_argument('--end-date',
        required=True,
        help=("End date in YYYY-MM-DD format"))
    parser.add_argument('--interval-days',
        default=7,
        type=int,
        help=("Number of days to stride by between --start-date and "
              "--end-date"))
    parser.add_argument('--cumulative',
        action='store_true',
        help=("If set, use all genomes up to each time point rather than "
              "only ones in the preceding window"))
    parser.add_argument('--sample-size-per-date-interval',
        default=100,
        type=int,
        help=("Number of sequences to sample randomly for each "
              "date, to increase in variation owing to having more "
              "sequences; if <= 0, do not subsample and use all sequences"))

    args = parser.parse_args()
    main(args)
