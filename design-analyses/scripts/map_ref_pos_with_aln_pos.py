"""Make map between position in alignment and a reference genome position.
"""

import argparse

from adapt import alignment
from adapt.utils import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


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


def main(args):
    seqs = seq_io.read_fasta(args.aln_fasta)

    # Find index of reference genome
    ref_idx = None
    seq_list = []
    for i, (seq_header, seq) in enumerate(seqs.items()):
        # Cut out accession in case it is [accession].[version]
        acc = seq_header.split(' ')[0].split('.')[0]
        if acc == args.ref_acc:
            assert ref_idx is None
            ref_idx = i
        seq_list += [seq]
    if ref_idx is None:
        raise Exception(("Reference accession was not found in FASTA"))

    # Make map of reference genome position to position in the alignment
    aln = alignment.Alignment.from_list_of_seqs(seq_list)
    m = ref_aln_pos_map(aln, ref_idx)

    # Write a map to TSV
    with open(args.out_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')
        header = ['ref-pos', 'aln-pos']
        write_row(header)
        for ref_pos in range(len(m)):
            row = [ref_pos, m[ref_pos]]
            write_row(row)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('aln_fasta',
        help=("Path to FASTA with alignment"))
    parser.add_argument('ref_acc',
        help=("Reference genome accession; must be a sequence "
              "header in the FASTA"))
    parser.add_argument('out_tsv',
        help=("TSV giving map between alignment and reference genome "
              "positions"))

    args = parser.parse_args()
    main(args)
