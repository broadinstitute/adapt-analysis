"""Script for randomly sampling sequences from a fasta file.
"""

import argparse
from collections import Counter
from collections import OrderedDict
import random

import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main(args):
    # Read input fasta
    seqs = seq_io.read_fasta(args.in_fasta)

    # Determine which sequence headers will be in the output
    names_to_sample = random.choices(list(seqs.keys()), k=args.num_seqs)

    # Count how often each sequence name should appear in the output
    name_count = Counter(names_to_sample)

    # A sequence name could have been sampled multiple times since
    # the sampling is with replacement, but the dicts are key'd on
    # sequence name (and the same sequence name should not show up
    # multiple times in the output fasta); replicate the multiplicity
    # here by renaming a sequence that should appear multiple times
    # from [name] to [name]-{1,2,...}
    seqs_sampled = OrderedDict()
    for name in sorted(set(names_to_sample)):
        if name_count[name] == 1:
            seqs_sampled[name] = seqs[name]
        else:
            # Add multiplicity to name
            for i in range(1, name_count[name] + 1):
                new_name = name + '-' + str(i)
                seqs_sampled[new_name] = seqs[name]

    # Write the resulting fasta
    seq_io.write_fasta(seqs_sampled, args.out_fasta)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('in_fasta')
    parser.add_argument('num_seqs', type=int,
            help=("Number of sequences to randomly sample with replacement"))
    parser.add_argument('out_fasta')

    args = parser.parse_args()
    main(args)
