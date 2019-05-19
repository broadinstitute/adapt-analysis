"""Read k-mers from sequence data and compute summary statistics.
"""

from collections import defaultdict
import glob
import itertools
import os

import seq_io

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def read_seqs(fasta_dir):
    """Read sequences from FASTAs in directory.

    Args:
        fasta_dir: path to directory containing .fasta.gz files

    Returns:
        dict {tax: sequences} where tax is a label for a taxonomy
        and sequences is an OrderedDict of sequences
    """
    seqs = {}
    for f in glob.glob(os.path.join(fasta_dir, '*.fasta.gz')):
        basename = os.path.basename(f)
        tax_name = basename.split('.')[0]
        seqs[tax_name] = seq_io.read_fasta(f)
    return seqs


def read_kmers(seqs, k=28):
    """Read k-mers from sequences

    Args:
        seqs: output of read_seqs()
        k: k-mer length

    Returns:
        tuple (x, y) where:
            x is a dict {k-mer: {(taxonomy identifier, sequence id)}} where
                the inner set gives pairs of (taxonomy, sequence) that
                contain the k-mer
            y is a dict {tax_name: taxonomy identifier} where taxonomy
                identifier is simply an integer
    """
    def kmers_from_seq(s):
        return [s[i:(i+k)] for i in range(len(s) - k + 1)]
    def is_unambig(kmer):
        return all(b in "ACGT" for b in kmer)

    seq_kmers = defaultdict(set)
    tax_ids = {}
    for tax_id, tax_name in enumerate(sorted(seqs.keys())):
        tax_ids[tax_name] = tax_id
        for seq_id, (_, seq) in enumerate(seqs[tax_name].items()):
            for kmer in kmers_from_seq(seq):
                seq_kmers[kmer].add((tax_id, seq_id))
    return (seq_kmers, tax_ids)


def compute_taxonomy_stats(seqs, out_tsv, k=28):
    """Determine statistics on each taxonomy.

    Args:
        seqs: output of read_seqs()
        out_tsv: path to TSV file to write values per taxonomy
        k: k-mer length
    """
    def kmers_from_seq(s):
        return [s[i:(i+k)] for i in range(len(s) - k + 1)]
    def is_unambig(kmer):
        return all(b in "ACGT" for b in kmer)

    # This stores a list of all k-mers to compute the number of unique
    # ones; there are much more efficient ways to do this (e.g., HLL)
    # but this only needs to be computed once
    rows = []
    all_kmers = []
    all_kmers_unambig = []
    for tax_name in sorted(seqs.keys()):
        num_seqs = len(seqs[tax_name])
        size = sum(len(seq) for _, seq in seqs[tax_name].items())
        kmers = list(itertools.chain.from_iterable(
                kmers_from_seq(seq) for _, seq in seqs[tax_name].items()))
        num_kmers = len(kmers)
        num_unique_kmers = len(set(kmers))
        all_kmers.extend(kmers)
        kmers_unambig = [kmer for kmer in kmers if is_unambig(kmer)]
        num_kmers_unambig = len(kmers_unambig)
        num_unique_kmers_unambig = len(set(kmers_unambig))
        all_kmers_unambig.extend(kmers_unambig)
        rows += [(tax_name, num_seqs, size, num_kmers, num_unique_kmers,
            num_kmers_unambig, num_unique_kmers_unambig)]

    print('Total number of k-mers:', len(all_kmers))
    print('Total number of unique k-mers:', len(set(all_kmers)))
    print('Total number of unambiguous k-mers:', len(all_kmers_unambig))
    print('Total number of unique unambiguous k-mers:',
            len(set(all_kmers_unambig)))

    with open(out_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join([str(x) for x in row]) + '\n')
        write_row(['tax_name', 'num_seqs', 'size', 'num_kmers',
            'num_unique_kmers', 'num_kmers_unambig', 'num_unique_kmers_unambig'])
        for row in rows:
            write_row(row)

