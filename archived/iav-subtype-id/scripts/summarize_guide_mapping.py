#!/usr/bin/env python

# Summarize mapping information of guides to target sequences.
#
# This reads in a SAM file of guides mapped to target sequences and
# outputs a summary of information about the fraction of target sequences
# that each guide (and guide set) map to.

import argparse
from collections import defaultdict
from collections import OrderedDict

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
            aln_info['query'] = ls[0]
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


def filter_aln_by_aln_score(aln, guides, score_thres,
        recompute_score=False):
    """Filter alignment entries by alignment score (AS).

    This also filters out alignments that are clipped, as determined
    based on the CIGAR string.

    Args:
        aln: output of read_sam()
        guides: dict giving guide sequence for each guide name
        score_thres: value of AS to keep
        recompute_score: if True, recompute an alignment score from
            the MD tag, in order to account for G-U pairing; otherwise,
            use the AS tag in the alignment

    Returns:
        aln with only entries whose alignment score value is >= score
    """
    aln_filt = defaultdict(list)
    for guide_name, aln_infos in aln.items():
        for aln_info in aln_infos:
            if 'S' in aln_info['cigar'] or 'H' in aln_info['cigar']:
                # This alignment is soft- or hard-clipped; filter it out by
                # skipping it
                continue

            if recompute_score:
                if 'MD' not in aln_info:
                    raise Exception("There is an alignment without MD field")
                if guide_name not in guides:
                    raise Exception(("There is a query with unknown guide "
                                     "sequence"))
                score = compute_aln_score(aln_info['MD'], guides[guide_name])
            else:
                if 'AS' not in aln_info:
                    raise Exception("There is an alignment without AS field")
                score = aln_info['AS']

            if score >= score_thres:
                aln_filt[guide_name].append(aln_info)
    return aln_filt


def read_guide_seqs(fn):
    """Read file listing guide sequence names.

    Args:
        fn: path to fasta file giving guide sequences (headers are
            used as names of guides)

    Returns:
        dict mapping guide names to guide sequences
    """
    guides = OrderedDict()
    with open(fn) as f:
        curr_guide_name = ""
        for line in f:
            line = line.rstrip()
            if len(line) == 0:
                # Reset the sequence being read on an empty line
                curr_guide_name = ""
                continue
            if curr_guide_name == "":
                # Must encounter a new sequence
                assert line.startswith('>')
            if line.startswith('>'):
                curr_guide_name = line[1:]
                assert curr_guide_name not in guides
                guides[curr_guide_name] = ""
            else:
                # Append the sequence
                guides[curr_guide_name] += line
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


def compute_aln_score(md_tag, read):
    """Compute an alignment score from the MD tag.

    This parses the MD tag of an alignment to calculate an alignment
    score. Here, we define the alignment score as the number of bases
    that "match" between the reference and guide sequence. We use
    "match" to mean if the bases are identical (in which case there will
    obviously be pairing between the guide and reference), as well as to
    account for G-U pairing (see next paragraph).

    An RNA guide with U can bind to target RNA with G, and vice-versa.
    Note that a guide sequence that ends up being synthesized to RNA
    is the reverse-complement of the guide sequence analyzed here, so
    that it will hybridize to the target (here, the guide sequences
    are designed to match the target). If the RNA guide were to have U and
    the RNA target were to have G, then the guide sequence here would be
    A and the reference would be G. If the RNA guide were to have G
    and the RNA target were to have U, then the guide sequence here were
    would be C and the reference would be T.

    Thus, we count a base X in the guide sequence as "matching" a
    base Y in the reference if any of the following are true:
      - X == Y
      - X == 'A' and Y == 'G'
      - X == 'C' and Y == 'T'

    We can count these using the MD tag and the guide sequence, based
    on the specs at:
      http://samtools.github.io/hts-specs/SAMtags.pdf
    and explained more at:
      https://github.com/vsbuffalo/devnotes/wiki/The-MD-Tag-in-BAM-Files

    To count the cases where X == Y, we just add up the integers in
    the tag. This works because:
      - mismatched bases are represented by the reference base and
        won't be counted
      - insertions are left out of the MD tag (they are not needed
        to determine the reference), so these won't be counted
      - deletions are represented by '^' followed by the reference
        base(s), so these won't be counted
      - soft clipped bases are left out of the MD tag, so they won't
        be counted
    (Note that, based on the arguments provided to bwa mem when
    making the alignments, the alignments should exclude insertions,
    deletions, and soft clipping.)

    To count the latter two cases, we look at the mismatched
    bases in the MD tag (for a mismatch between the read (i.e., guide)
    and reference, the reference base is in the MD tag) and compare it
    at that position to the base in the read. We can do this easily if
    we assume (and verify) that there are no indels in the alignment
    between the guide and reference.

    Note that this will fail (namely, an assertion at the end) if
    the alignment is soft- or hard-clipped.

    Args:
        md_tag: MD tag for an alignment as given in a SAM file
        read: read sequence in the corresponding alignment (in this
            case, the guide sequence)

    Returns:
        integer representing an alignment score, calculated as
        described above
    """
    curr_pos = 0
    curr_matching_length = ''
    num_matching = 0

    def process_end_of_matching_length():
        nonlocal curr_pos
        nonlocal curr_matching_length
        nonlocal num_matching
        # curr_matching_length is a string giving a number of
        # matching bases
        matching_length = int(curr_matching_length)
        curr_pos += matching_length
        curr_matching_length = ''
        num_matching += matching_length

    for c in md_tag:
        if '0' <= c <= '9':
            # Indicates part of an integer giving a number of matching bases
            curr_matching_length += c
        else:
            # Check the character is not '^', which indicates a deletion
            if c == '^':
                raise Exception("MD tag indicates deletion")
            if len(curr_matching_length) > 0:
                process_end_of_matching_length()

            # The character should be A, T, C, or G and indicates a
            # mismatch between the guide (read) and reference
            assert c in ['A', 'T', 'C', 'G']
            ref_base = c
            guide_base = read[curr_pos]

            # Count a match if there is G-U pairing
            if guide_base == 'A' and ref_base == 'G':
                num_matching += 1
            if guide_base == 'C' and ref_base == 'T':
                num_matching += 1

            # Advance the current position (used up by this mismatch)
            curr_pos += 1
    if len(curr_matching_length) > 0:
        process_end_of_matching_length()

    # Verify that we went through all bases in the read (i.e., there are
    # no deletions or insertions relative to the reference; without this
    # check we might not be able to tell if there was an insertion, in
    # which case the G-U pairing checks above would have been incorrect)
    if curr_pos != len(read):
        print(md_tag, curr_pos, len(read), read)
    assert curr_pos == len(read)

    return num_matching


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
    guides = read_guide_seqs(args.guide_seqs)

    if args.aln_score_filter:
        aln = filter_aln_by_aln_score(aln, guides, args.aln_score_filter,
            recompute_score=args.recompute_score_for_filter)

    # Find the set of sequences that each guide hits
    seqs_hit = seqs_mapped_per_guide(aln)

    # Calculation fraction of target sequences hit by each guide and
    # guide set
    frac_hit_per_guide = frac_of_refs_hit_by_guide(seqs_hit,
        guides.keys(), args.num_refs)
    frac_hit_per_guide_set = frac_of_refs_hit_by_guide_set(seqs_hit,
        guides.keys(), args.num_refs)

    # Write results to TSV
    write_frac_hits(frac_hit_per_guide, args.out_for_guide)
    write_frac_hits(frac_hit_per_guide_set, args.out_for_guide_set)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('in_sam', help='Path to input SAM')
    parser.add_argument('guide_seqs',
        help=('Path to fasta giving all guide sequences (header names giving '
              'guide names'))
    parser.add_argument('num_refs', type=int, help='Number of target sequences')
    parser.add_argument('out_for_guide',
        help='Path to output TSV with info for each guide')
    parser.add_argument('out_for_guide_set',
        help='Path to output TSV with info for each guide set')

    parser.add_argument('--aln-score-filter', type=int,
        help=('Filter on this alignment score (>= value)'))
    parser.add_argument('--recompute-score-for-filter', action='store_true',
        help=('When computing alignment score for filtering, recompute '
              'it using the MD field in order to account for G-U pairing; '
              'otherwise, use the AS field value in the alignment'))

    args = parser.parse_args()
    main(args)
