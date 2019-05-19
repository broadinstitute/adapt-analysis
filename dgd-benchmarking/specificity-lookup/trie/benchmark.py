"""Benchmark trie data structure.
"""

import read_kmers

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def main():
    seqs = read_kmers.read_seqs('../data/fastas')
    read_kmers.compute_taxonomy_stats(seqs, 'out/stats.tsv')

    kmers, tax_ids = read_kmers.read_kmers(seqs)
    print(tax_ids)
    print(kmers)


if __name__ == "__main__":
    main()
