"""Functions for computing statistics from one or more designs.
"""

import argparse
from collections import defaultdict
import glob
import math
import statistics

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
        Shannon entropy
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
        entropy += -1.0 * frac * math.log2(frac)

    return entropy


def compute_jensen_shannon_divergence(designs1, designs2):
    """Compute Jensen-Shannon divergence between two collections of designs.

    One option is to compute KL divergence between P(x) and Q(x) where P(x)
    represents a probability distribution over the designs in designs1,
    and likewise for Q(x) for designs2. However, this would be difficult here
    (without smoothing the distributions) because the probability distributions
    contain many zeros: there are values x (designs) where Q(x) = 0 but
    P(x) != 0. Due to these values, the KL divergence is undefined.

    Jensen-Shannon divergence is appealing because it does not have the
    above problem, it is symmetric (KL divergence is asymmetric), and simply
    taking the square root turns it into a distance metric. Also, if we
    use base 2 for the logarithm, its value is in [0,1].

    We can compute the Jensen-Shannon divergence (JSD) as follows:
    Let P represent a probability distribution for designs1, Q represent
    a probability distribution for designs2, and M be (1/2)*(P + Q).
    By definition,
        JSD = (1/2)*sum(P(x)*log(P(x)/M(x))) + (1/2)*sum(Q(x)*log(Q(x)/M(x)))
    where sums are taken over designs x for x in the union of designs1
    and designs2. Simplifying this:
        JSD = (1/2)*(-H(P) - H(Q) - sum((P(x) + Q(x))*log(M(x))))
            = -(1/2)*sum((P(x) + Q(x))*log(M(x))) - (1/2)*H(P) - (1/2)*H(Q)
            = -sum(M(x)*log(M(x))) - (1/2)*(H(P) + H(Q))
            = H(M) - (1/2)*(H(P) + H(Q))
    where H(.) is the Shannon entropy.

    Args:
        designs1: collection of Design objects
        designs2: collection of Design objects

    Returns:
        Jensen-Shannon divergence
    """
    # Compute the entropy of M = (1/2)*(P + Q)
    P_design_count = defaultdict(int)
    for design in designs1:
        P_design_count[design] += 1
    P_design_frac = {design: float(count) / len(designs1)
            for design, count in P_design_count.items()}
    Q_design_count = defaultdict(int)
    for design in designs2:
        Q_design_count[design] += 1
    Q_design_frac = {design: float(count) / len(designs2)
            for design, count in Q_design_count.items()}
    M_entropy = 0.0
    for design in set(designs1 + designs2):
        if design in P_design_frac:
            P_prob = P_design_frac[design]
        else:
            P_prob = 0
        if design in Q_design_frac:
            Q_prob = Q_design_frac[design]
        else:
            Q_prob = 0
        M_prob = 0.5 * (P_prob + Q_prob)
        M_entropy += -1.0 * M_prob * math.log2(M_prob)

    # Compute the entropy of P and Q
    P_entropy = compute_entropy(designs1)
    Q_entropy = compute_entropy(designs2)

    jsd = M_entropy - 0.5 * (P_entropy + Q_entropy)
    return jsd


def compute_pairwise_jaccard_similarity_within_designs(designs):
    """Compute pairwise Jaccard similarities of a collection of designs.

    Args:
        designs: collection of Design objects

    Returns:
        list of pairwise Jaccard similarity between all pairs of designs
        within a collection
    """
    if len(designs) == 1:
        # Undefined
        return None

    similarities = []
    for i in range(len(designs)):
        for j in range(i+1, len(designs)):
            similarities += [designs[i].jaccard_similarity(designs[j])]
    return similarities


def compute_pairwise_jaccard_similarity_between_designs(designs1,
        designs2):
    """Compute pairwise Jaccard similarity between two collections
    of designs.

    Note that if designs1 == designs2, the value returned by this function
    should be slightly greater than the value returned by
    compute_average_pairwise_jaccard_similarity_within_designs(designs1),
    even if it seems that they should be equal. The reason is that this
    function will include, in computing the average, pairs of the
    same design (i.e., designs1[i] and designs1[j] where i==j), whereas
    the function averaging *within* designs does *not* include, in computing
    its average, identical designs (for which the Jaccard similarity is
    1.0).

    Args:
        designs1: collection of Design objects
        designs2: collection of Design objects

    Returns:
        list of pairwise Jaccard similarity between all pairs of designs
        across two collections
    """
    similarities = []
    for i in range(len(designs1)):
        for j in range(len(designs2)):
            similarities += [designs1[i].jaccard_similarity(designs2[j])]
    return similarities


def run_dispersion(args):
    designs = read_designs(args.design_tsvs, num_targets=args.num_targets)

    entropy = compute_entropy(designs)
    jaccard_similarities = compute_pairwise_jaccard_similarity_within_designs(designs)

    if args.pairwise_jaccard_distribution_out:
        with open(args.pairwise_jaccard_distribution_out, 'w') as f:
            for s in jaccard_similarities:
                f.write(str(s) + '\n')

    print("Entropy:", entropy)
    print("Average pairwise Jaccard similarity within collection of designs:",
            statistics.mean(jaccard_similarities))


def run_compare(args):
    designs1 = read_designs(args.design_tsvs1, num_targets=args.num_targets)
    designs2 = read_designs(args.design_tsvs2, num_targets=args.num_targets)

    jsd = compute_jensen_shannon_divergence(designs1, designs2)
    jaccard_similarities = compute_pairwise_jaccard_similarity_between_designs(
            designs1, designs2)

    if args.pairwise_jaccard_distribution_out:
        with open(args.pairwise_jaccard_distribution_out, 'w') as f:
            for s in jaccard_similarities:
                f.write(str(s) + '\n')

    print("Jensen-Shannon divergence:", jsd)
    print("Average pairwise Jaccard similarity between collections of designs:",
            statistics.mean(jaccard_similarities))


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
    parser_dispersion.add_argument('--pairwise-jaccard-distribution-out',
            help=("If specified, path to which to write the distribution "
                  "of pairwise Jaccard similarity values"))

    parser_compare = subparsers.add_parser('compare',
            parents=[common_parser],
            help=("Compute stats comparing two collections of designs."))
    parser_compare.add_argument('--design-tsvs1',
            nargs='+',
            help=("Paths to one of more TSV files containing designs; use "
                  "* or ** as wildcards"))
    parser_compare.add_argument('--design-tsvs2',
            nargs='+',
            help=("Paths to one of more TSV files containing designs; use "
                  "* or ** as wildcards"))
    parser_compare.add_argument('--pairwise-jaccard-distribution-out',
            help=("If specified, path to which to write the distribution "
                  "of pairwise Jaccard similarity values"))

    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if args.command == 'dispersion':
        run_dispersion(args)
    elif args.command =='compare':
        run_compare(args)
    else:
        raise Exception(("Unknown command %s") % args.command)


if __name__ == "__main__":
    main()
