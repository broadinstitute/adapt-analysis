#!/usr/bin/env python

# Summarize mapping information of guides to target sequences.
#
# This reads in a SAM file of guides mapped to target sequences and
# outputs a summary of information about the fraction of target sequences
# that each guide (and guide set) map to.

import argparse
from collections import defaultdict

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def read_sam(in_sam):
    """Map each guide name to a list of dicts describing the alignment.

    Args:
        in_sam: path to input SAM file

    Returns:
        dict {guide name: [aligments]}, value is list in case a guide
        aligns more than once to the target sequences
    """
    aln = defaultdict(list)
    with open(in_sam) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('@'):
                # Skip header
                continue
            ls = line.split('\t')
            guide_name = ls[0]

            aln_info = {}
            aln_info['ref'] = ls[2]
            aln_info['pos'] = int(ls[3])
            aln_info['cigar'] = ls[5]

            # Read optional fields
            opt_fields = {}
            for opt_field in ls[11:]:            
                k, t, v = opt_field.split(':')
                if t == 'i':
                    # Value of field is an int
                    v = int(v)
                opt_fields[k] = v
            for k in ('NM', 'AS', 'MD'):
                if k in opt_fields:
                    aln_info[k] = opt_fields[k]

            aln[guide_name].append(aln_info)
    return aln


def filter_aln_by_aln_score(aln, score):
    """Filter alignment entries by alignment score (AS).

    Args:
        aln: output of read_sam()
        score: value of AS to keep

    Returns:
        aln with only entries whose AS value is >= score
    """
    aln_filt = defaultdict(list)
    for guide_name, aln_infos in aln.items():
        for aln_info in aln_infos:
            if 'AS' not in aln_info:
                raise Exception("There is an alignment without AS field")
            if aln_info['AS'] >= score:
                aln_filt[guide_name].append(aln_info)
    return aln_filt


def read_list_of_guide_names(guide_list):
    """Read file listing guide sequence names.

    Args:
        guide_list: path to file in which each line lists a guide header

    Returns:
        list of lines in the file
    """
    guides = []
    with open(guide_list) as f:
        for line in f:
            guides += [line.rstrip()]
    return guides


def seqs_mapped_per_guide(aln):
    """Find set of ref sequences mapped to for each guide.

    Args:
        aln: output of read_sam()

    Returns:
        dict mapping guide name to set of sequences that guide aligns to
    """
    hits = defaultdict(set)
    for guide_name, aln_infos in aln.items():
        for aln_info in aln_infos:
            hits[guide_name].add(aln_info['ref'])
    return hits


def frac_of_refs_hit_by_guide(hits, guides, num_refs):
    """Calculate fraction of reference (target) sequences each guide aligns to.

    Args:
        hits: dict mapping guide name to set of sequences that the guide
            aligns to; this should be a defaultdict so that, if a guide
            is not present, it yields an empty set
        guides: list of all guide sequence names
        num_refs: number of reference (target) sequences

    Returns:
        dict mapping each guide in guides to the fraction of sequences
        (denominator is num_refs) that the guide hits
    """
    frac_hit = {}
    for guide in guides:
        num_hit = len(hits[guide])
        frac_hit[guide] = float(num_hit) / num_refs
    return frac_hit


def frac_of_refs_hit_by_guide_set(hits, guides, num_refs):
    """Calculate fraction of reference (target) sequences each guide set
    aligns to.

    Guide sets are determined from the guide names. A guide name should
    be in the format:
        [subtype]-u[N]-g[M]
    where each u[N] defines a set of guides.

    A guide set hits a sequence if any of the guides in it aligns to the
    sequence.

    Args:
        hits: dict mapping guide name to set of sequences that the guide
            aligns to; this should be a defaultdict so that, if a guide
            is not present, it yields an empty set
        guides: list of all guide sequence names
        num_refs: number of reference (target) sequences

    Returns:
        dict mapping each guide set (determined from guide names in guides)
        to the fraction of sequences (denominator is num_refs) that the
        the set of guides collectively hits
    """
    # Find guide sets and the guides in each
    guide_sets = defaultdict(set)
    for guide in guides:
        guide_split = guide.split('-')
        assert len(guide_split) == 3

        assert guide_split[1][0] == 'u'
        assert guide_split[2][0] == 'g'

        guide_set = guide_split[0] + '-' + guide_split[1]
        guide_sets[guide_set].add(guide)

    frac_hit = {}
    for guide_set in guide_sets.keys():
        # Find the sequences collectively covered by the guides in guide_set
        all_hit = set()
        for guide in guide_sets[guide_set]:
            all_hit.update(hits[guide])
        num_hit = len(all_hit)
        frac_hit[guide_set] = float(num_hit) / num_refs
    return frac_hit


def write_frac_hits(frac_hit, out_tsv):
    """Write dict giving fraction of sequences hit to a TSV.

    Args:
        frac_hit: dict where key is a guide or guide set and the value
            is a fraction of target sequences hit by the key
        out_tsv: path to output TSV file
    """
    with open(out_tsv, 'w') as f:
        for k in sorted(frac_hit.keys()):
            f.write('\t'.join([k, str(frac_hit[k])]) + '\n')


def main(args):
    aln = read_sam(args.in_sam)
    guides = read_list_of_guide_names(args.guide_list)

    if args.aln_score_filter:
        aln = filter_aln_by_aln_score(aln, args.aln_score_filter)

    # Find the set of sequences that each guide hits
    seqs_hit = seqs_mapped_per_guide(aln)

    # Calculation fraction of target sequences hit by each guide and
    # guide set
    frac_hit_per_guide = frac_of_refs_hit_by_guide(seqs_hit,
        guides, args.num_refs)
    frac_hit_per_guide_set = frac_of_refs_hit_by_guide_set(seqs_hit,
        guides, args.num_refs)

    # Write results to TSV
    write_frac_hits(frac_hit_per_guide, args.out_for_guide)
    write_frac_hits(frac_hit_per_guide_set, args.out_for_guide_set)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('in_sam', help='Path to input SAM')
    parser.add_argument('guide_list', help='Path to list of all guides')
    parser.add_argument('num_refs', type=int, help='Number of target sequences')
    parser.add_argument('out_for_guide',
        help='Path to output TSV with info for each guide')
    parser.add_argument('out_for_guide_set',
        help='Path to output TSV with info for each guide set')

    parser.add_argument('--aln-score-filter', type=int,
        help=('Filter on this alignment score (>= value)'))

    args = parser.parse_args()
    main(args)
