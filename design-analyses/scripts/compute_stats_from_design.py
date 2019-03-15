"""Functions for computing statistics from one or more designs.
"""

import argparse
from collections import defaultdict
import glob
import math

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class DesignTarget:
    """Store information on a design of a single target.

    To enable a comparison of designs across different alignments, this
    does not store positional information (only sequence).
    """

    def __init__(self, guide_seqs, left_primer_seqs, right_primer_seqs,
            cost):
        self.guide_seqs = tuple(sorted(guide_seqs))
        self.left_primer_seqs = tuple(sorted(left_primer_seqs))
        self.right_primer_seqs = tuple(sorted(right_primer_seqs))
        self.cost = cost

    def __eq__(self, other):
        return (self.guide_seqs == other.guide_seqs and
                self.left_primer_seqs == other.left_primer_seqs and
                self.right_primer_seqs == other.right_primer_seqs)

    def __hash__(self):
        return (hash(self.guide_seqs) +
                hash(self.left_primer_seqs) +
                hash(self.right_primer_seqs))


class Design:
    """Store information on a design encompassing multiple possible targets.

    This stores the targets as a set (unordered). As a result, two designs
    can be equal if all their targets are equal, even if they are ordered
    differently (e.g., different costs for each).
    """

    def __init__(self, targets):
        self.targets = frozenset(set(targets))

    def __eq__(self, other):
        return self.targets == other.targets

    def __hash__(self):
        return hash(self.targets)

    def jaccard_similarity(self, other):
        """Compute Jaccard similarity between this design and another.

        Args:
            other: Design object

        Returns:
            Jaccard similarity, by comparing targets, between self and other
        """
        intersection = self.targets & other.targets
        union = self.targets | other.targets
        return float(len(intersection)) / float(len(union))

    @staticmethod
    def from_file(fn, num_targets=None):
        """Read a collection of targets from a file.

        Args:
            fn: path to a TSV file giving targets
            num_targets: only construct a Design from the top num_targets
                targets, as ordered by cost (if None, use all)

        Returns:
            Design object
        """
        rows = []
        with open(fn) as f:
            col_names = {}
            for i, line in enumerate(f):
                line = line.rstrip()
                ls = line.split('\t')
                if i == 0:
                    # Parse header
                    for j in range(len(ls)):
                        col_names[j] = ls[j]
                else:
                    # Read each column as a variable
                    cols = {}
                    for j in range(len(ls)):
                        cols[col_names[j]] = ls[j]
                    rows += [(cols['cost'], cols['target-start'],
                             cols['target-end'], cols)]

        # Sort rows by cost (first in the tuple); in case of ties, sort
        # by target start and target end positions (second and third in
        # the tuple)
        # Pull out the best N targets
        rows = sorted(rows)
        if num_targets != None:
            if len(rows) < num_targets:
                raise Exception(("The number of rows in a design (%d) is fewer "
                    "than the number of targets to read (%d)") %
                    (len(rows), num_targets))
            rows = rows[:num_targets]

        targets = []
        for row in rows:
            _, _, _, cols = row
            targets += [DesignTarget(
                cols['guide-target-sequences'].split(' '),
                cols['left-primer-target-sequences'].split(' '),
                cols['right-primer-target-sequences'].split(' '),
                float(cols['cost'])
            )]

        return Design(targets)


def read_designs(design_tsvs, num_targets=None):
    """Read a collection of designs.

    Args:
        design_tsvs: paths to one or more TSV files containing designs;
            can contain * or ** as wildcards
        num_targets: only construct a Design from the top num_targets
            targets, as ordered by cost (if None, use all)

    Returns:
        collection of Design objects
    """
    designs = []
    for design_tsv in design_tsvs:
        # Use glob to expand wildcards
        files = glob.glob(design_tsv)
        for fn in files:
            designs += [Design.from_file(fn, num_targets=num_targets)]
    return designs


def compute_entropy(designs):
    """Compute Shannon entropy of a collection of designs.

    Args:
        designs: collection of Design objects

    Returns:
        Shannon entroy
    """
    # Count the unique designs
    design_count = defaultdict(int)
    for design in designs:
        design_count[design] += 1

    # Compute entropy by summing over the fraction (probability) of
    # each unique design
    num_designs = len(designs)
    entropy = 0.0
    for design, count in design_count.items():
        frac = float(count) / num_designs
        entropy += -1.0 * frac * math.log(frac)

    return entropy


def compute_average_pairwise_jaccard_similarity(designs):
    """Compute average pairwise Jaccard similarity of designs.

    Args:
        designs: collection of Design objects

    Returns:
        average pairwise Jaccard similarity
    """
    if len(designs) == 1:
        return 1.0

    total_similarity = 0
    n = 0
    for i in range(len(designs)):
        for j in range(i, len(designs)):
            total_similarity += designs[i].jaccard_similarity(designs[j])
            n += 1
    return float(total_similarity) / n


def run_dispersion(args):
    designs = read_designs(args.design_tsvs, num_targets=args.num_targets)

    entropy = compute_entropy(designs)
    jaccard_similarity = compute_average_pairwise_jaccard_similarity(designs)

    print("Entropy:", entropy)
    print("Average pairwise Jaccard similarity:", jaccard_similarity)


def parse_args():
    parser = argparse.ArgumentParser()

    common_parser = argparse.ArgumentParser(add_help=False)
    common_parser.add_argument('--num-targets',
            type=int,
            default=10,
            help=("Only read the top NUM_TARGETS targets (according to cost)"))

    subparsers = parser.add_subparsers(
            title='commands',
            dest='command')

    parser_dispersion = subparsers.add_parser('dispersion',
            parents=[common_parser],
            help=("Compute dispersion of a collection of designs."))
    parser_dispersion.add_argument('design_tsvs',
            nargs='+',
            help=("Paths to one of more TSV files containing designs; use "
                  "* or ** as wildcards"))

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if args.command == 'dispersion':
        run_dispersion(args)
    else:
        raise Exception(("Unknown command %s") % args.command)


if __name__ == "__main__":
    main()
