"""Select species from NCBI's viral genome neighbor accession list.
"""

import argparse
from collections import defaultdict
import gzip

__author__ = 'Hayden Metsky <hayden@mit.edu>'


class SequenceFromAccessionList:
    def __init__(self, representative, accession, host,
                 lineage, taxonomy_name, segment):
        self.representative = representative
        self.accession = accession
        self.hosts = host.split(',')
        self.lineage = lineage
        self.taxonomy_name = taxonomy_name
        self.segment = segment
        self.human_is_host = 'human' in host
        self.vertebrate_is_host = ('vertebrate' in self.hosts or 'human' in self.hosts
                or 'vertebrates' in self.hosts or 'humans' in self.hosts)
        self.is_segmented = segment != 'segment'

        self.segment = self.segment.replace('segment ', '').replace('segment', '')
        if self.segment == '':
            self.segment = 'None'

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    @staticmethod
    def from_line(line):
        ls = line.split('\t')
        if ls[2] == '':
            # In a recent accession list, for some accessions host is blank;
            # ignore these
            return None
        return SequenceFromAccessionList(ls[0], ls[1], ls[2], ls[3], ls[4], ls[5])


def read_genome_accession_list(fn, force_exclude_taxa=None):
    """Read genome accession list.

    Args:
        fn: path to accession list

    Returns:
        list of SequenceFromAccessionList objects
    """
    # Avoid reading an accession twice; it may appear twice because of
    # merged lists
    read_accessions = set()

    if force_exclude_taxa is None:
        force_exclude_taxa = set()

    sequences = []
    with gzip.open(fn, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            s = SequenceFromAccessionList.from_line(line.rstrip())
            if s is None:
                continue
            if s.lineage in force_exclude_taxa:
                continue
            if s.accession in read_accessions:
                continue
            else:
                sequences += [s]
                read_accessions.add(s.accession)
    return sequences


def read_force_include_taxa(fn):
    """Read taxa to include regardless of how host is labeled.

    Args:
        fn: path to tsv

    Returns:
        set of lineages
    """
    lineages = set()
    with open(fn) as f:
        for i, line in enumerate(f):
            if i == 0:
                # Header
                continue
            ls = line.rstrip().split('\t')
            lineage = ','.join([ls[0], ls[1], ls[2]])
            lineages.add(lineage)
    return lineages


def read_force_exclude_taxa(fn):
    """Read taxa to exclude

    Args:
        fn: path to tsv

    Returns:
        set of lineages
    """
    lineages = set()
    with open(fn) as f:
        for i, line in enumerate(f):
            if i == 0:
                # Header
                continue
            ls = line.rstrip().split('\t')
            lineage = ','.join([ls[0], ls[1], ls[2]])
            lineages.add(lineage)
    return lineages


def lineage_refseqs(sequences):
    """Find RefSeqs for each lineage.

    Args:
        sequences: list of SequenceFromAccessionList objects

    Returns:
        dict {lineage: set of RefSeq}
    """
    r = defaultdict(set)
    for s in sequences:
        r[s.lineage].add(s.representative)
    return r


def filter_to_have_family(sequences):
    """Filter sequences to only contain ones with a family.

    Args:
        sequences: list of SequenceFromAccessionList objects
    
    Returns:
        list of SequenceFromAccessionList objects
    """
    def has_family(lineage):
        if ',' not in lineage:
            return False
        ls = lineage.split(',')
        if len(ls) == 1:
            return False
        if len(ls[0]) == 0:
            return False
        return True

    return [s for s in sequences if has_family(s.lineage)]


def count_seqs_per_lineage(sequences):
    """Count number of sequences for each segment for each lineage.

    Args:
        sequences: list of SequenceFromAccessionList objects

    Returns:
        dict {lineage: {segment: number of sequences}}
    """
    c = {}
    for s in sequences:
        if s.lineage not in c:
            c[s.lineage] = defaultdict(int)
        c[s.lineage][s.segment] += 1
    return c


def filter_by_host(sequences, force_include_taxa=None, force_exclude_taxa=None):
    """Filter sequences to only ones with vertebrate host.

    Args:
        sequences: list of SequenceFromAccessionList objects
        force_include_taxa: if set, set of lineages to include regardless of
            how host is labeled
        force_exclude_taxa: if set, set of lineages to exclude
    
    Returns:
        list of SequenceFromAccessionList objects
    """
    if force_include_taxa is None:
        force_include_taxa = set()
    if force_exclude_taxa is None:
        force_exclude_taxa = set()
    return [s for s in sequences if
            (s.vertebrate_is_host or s.lineage in force_include_taxa) and s.lineage not in force_exclude_taxa]


def pick_segment(sequences):
    """Pick a segment for each lineage.

    Use the one with the most sequences, except for IAV and IBV. The latter
    two are outliers: use segment 2 for IAV and segment 1 for IBV (the
    most conserved).

    The reason for picking the segment with the most sequences is just
    that these may be the segments for which we know the most about
    sequence diversity. On the other hand, they may also be the ones
    with the most sequence diversity and thus the hardest to detect.

    Args:
        sequences: list of SequenceFromAccessionList objects
    
    Returns:
        list of SequenceFromAccessionList objects
    """
    counts = count_seqs_per_lineage(sequences)

    segment_chosen = {}
    for lineage in counts.keys():
        if 'Influenza A virus' in lineage:
            segment_chosen[lineage] = '2'
        elif 'Influenza B virus' in lineage:
            segment_chosen[lineage] = '1'
        else:
            most_common = sorted(list(counts[lineage].items()),
                    key=lambda x: x[1], reverse=True)[0][0]
            segment_chosen[lineage] = most_common
    return segment_chosen


def main(args):
    if args.force_include_taxa:
        forced_include_taxa = read_force_include_taxa(args.force_include_taxa)
    else:
        forced_include_taxa = set()
    if args.force_exclude_taxa:
        forced_exclude_taxa = read_force_exclude_taxa(args.force_exclude_taxa)
    else:
        forced_exclude_taxa = set()

    sequences = read_genome_accession_list(args.accession_list,
            force_exclude_taxa=forced_exclude_taxa)

    sequences = filter_by_host(sequences, force_include_taxa=forced_include_taxa,
            force_exclude_taxa=forced_exclude_taxa)
    sequences = filter_to_have_family(sequences)

    segments = pick_segment(sequences)
    refseqs = lineage_refseqs(sequences)
    counts = count_seqs_per_lineage(sequences)

    lineages = sorted(list(set([s.lineage for s in sequences])))
    with open(args.out_tsv, 'w') as fw:
        header = ['family', 'genus', 'species', 'segment', 'refseqs', 'neighbor-count']
        fw.write('\t'.join(header) + '\n')

        for lineage in lineages:
            ls = lineage.split(',')
            if len(ls) == 2:
                # For simplicity, say this is family and species (even though
                # it may actually be genus, not family)
                family, species = ls
                genus = ''
            else:
                family, genus, species = ls
            if genus == '':
                genus = 'unknown'
            rs = ','.join(list(refseqs[lineage]))
            s = segments[lineage]
            neighbor_count = counts[lineage][s]
            if 'Influenza A virus' in species or 'Influenza B virus' in species:
                # The accession list does not store their accessions
                neighbor_count = 'NA'

            row = [family, genus, species, s, rs, neighbor_count]

            fw.write('\t'.join(str(x) for x in row) + '\n')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('accession_list',
        help=("Path to genome accession list, gzip'd"))
    parser.add_argument('out_tsv',
        help=("Path to output TSV"))

    parser.add_argument('--force-include-taxa',
        help=("Path to TSV of taxa to include regardless of how host is labeled"))
    parser.add_argument('--force-exclude-taxa',
        help=("Path to TSV of taxa to exclude because they are duplicates of another"))

    args = parser.parse_args()
    main(args)
