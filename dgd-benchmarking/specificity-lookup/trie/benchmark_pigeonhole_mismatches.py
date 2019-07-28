"""Benchmark trie data structure, with queries based on ensuring there are
floor(m/2) mismatches in each half of the 28-mer.

The basic idea here is: we store two tries. For each 28-mer x, we place a
28-mer x into one trie and the reverse of x into the other trie. if we want
to query a 28-mer q up to m mismatches, then by the pigeonhole principle
at least one of the two halves must have at most floor (m/2) mismatches. So
we can query q in the first try, tolerating only floor(m/2) mismatches for the
first 14 levels. And we can query the reverse of q in the other trie,
tolerating only floor(m/2) mismatches for the first 14 levels.
"""

from collections import defaultdict
import logging
import math
import random
import time

import benchmark
import read_kmers
import trie

__author__ = 'Hayden Metsky <hayden@mit.edu>'

# Configure basic logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def query_for_taxonomy(t_0, t_1, taxid, kmer_sample_size=100):
    """Query random sampling of k-mers from a given taxonomy in the trie
    spaces.

    This first masks that taxonomy from the tries (otherwise most results
    will be for that taxonomy).

    Args:
        t_0: trie.Trie object with forward k-mers
        t_1: trie.Trie object with reverse k-mers
        taxid: taxonomy identifier
        kmer_sample_size: number of k-mers to sample for querying

    Returns:
        (num_matches, num_nodes_visited, runtime) where
        each is a dict {(gu_pairing, mismatches): v} where v is a list of
        values across the randomly sampled k-mers
    """
    # Find all the k-mers in taxid
    logging.info("Finding k-mers for taxid %d", taxid)
    leaves = t_0.root_node.traverse_and_find(taxid)

    # Construct a list of k-mers and a corresponding list of the number of
    # sequences in taxid that contain each k-mer, to use for weighting
    # during sampling
    taxid_kmers = []
    taxid_kmer_occ = []
    for kmer, li in leaves:
        num_seqs_with_kmer_in_taxid = len(li.d[taxid])
        taxid_kmers += [kmer]
        taxid_kmer_occ += [num_seqs_with_kmer_in_taxid]

    if len(taxid_kmers) == 0:
        logging.info("No k-mers for taxid %d", taxid)
        return None

    # Randomly sample k-mers with replacement, weighting the selection by
    # the number of sequences in taxid that contain each k-mer
    taxid_kmers_sample = random.choices(taxid_kmers,
            weights=taxid_kmer_occ,
            k=kmer_sample_size)

    # Mask taxid from the tries
    logging.info("Masking taxid %d", taxid)
    t_0.mask(taxid)
    t_1.mask(taxid)

    logging.info("Querying for taxid %d", taxid)

    # Query each k-mer at different parameters (but all with G-U pairing)
    num_matches = {}
    perf_num_nodes_visited = {}
    perf_runtime = {}
    for gu_pairing in [False, True]:
        for m in [0, 1, 2, 3]:
            num_matches_for_params = []
            perf_num_nodes_visited_for_params = []
            perf_runtime_for_params = []
            for kmer in taxid_kmers_sample:
                start_time = time.time()

                # By the pigeonhole principle, we'll have at most floor(m/2)
                # mismatches in one of the two halves of the 28-mer
                m_half = int(m / 2)
                assert m_half in [0, 1]

                num_results = 0
                num_nodes_visited = 0

                # Query for tries maybe containing each half of the kmer
                # Query the forward and reverse in the corresponding trie
                for reverse, t in [(False, t_0), (True, t_1)]:
                    if reverse:
                        # Query kmer in reverse
                        q = kmer[::-1]
                    else:
                        q = kmer
                    # Only allow up to m_half mismatches for the first
                    # half of the search (i.e., levels 0 through 13)
                    q_results, q_num_nodes_visited = t.query(
                            q, mismatches=m, gu_pairing=gu_pairing,
                            mismatches_to_level=(m_half, 13))
                    num_results += len(q_results)
                    num_nodes_visited += q_num_nodes_visited

                end_time = time.time()

                num_matches_for_params += [num_results]
                perf_num_nodes_visited_for_params += [num_nodes_visited]
                perf_runtime_for_params += [end_time - start_time]

            num_matches[(gu_pairing, m)] = num_matches_for_params
            perf_num_nodes_visited[(gu_pairing, m)] = perf_num_nodes_visited_for_params
            perf_runtime[(gu_pairing, m)] = perf_runtime_for_params

    # Unmask taxid
    logging.info("Unmasking taxid %d", taxid)
    t_0.unmask_all()
    t_1.unmask_all()

    return (num_matches,
            perf_num_nodes_visited, perf_runtime)


def benchmark_queries_across_taxonomies(t_0, t_1, tax_ids, out_tsv,
        tax_ids_to_sample=None):
    """Benchmark queries across taxonomies.

    Args:
        t_0: trie.Trie object with forward k-mers
        t_1: trie.Trie object with reverse k-mers
        tax_ids: dict {tax_name: taxonomy identifier} where taxonomy
            identifier is simply an integer
        out_tsv: path to TSV file to which to write benchmark results
        tax_ids_to_sample: number of tax ids to sample
    """
    if tax_ids_to_sample is not None:
        tax_ids_to_insert = set(random.sample(tax_ids.keys(), tax_ids_to_sample))
        tax_ids_subsampled = {}
        for k in tax_ids_to_insert:
            tax_ids_subsampled[k] = tax_ids[k]
        tax_ids = tax_ids_subsampled

    benchmark_results = []
    for tax_name, taxid in tax_ids.items():
        x = query_for_taxonomy(t_0, t_1, taxid)
        if x is None:
            # No k-mers for taxid
            continue
        num_matches, perf_num_nodes_visited, perf_runtime = x
        for benchmark_name, d in [('matches', num_matches),
                ('nodes_visited', perf_num_nodes_visited),
                ('runtime', perf_runtime)]:
            for gu_pairing, mismatches in d.keys():
                for v in d[(gu_pairing, mismatches)]:
                    benchmark_results += [(tax_name, benchmark_name,
                        gu_pairing, mismatches, v)]

    with open(out_tsv, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join([str(x) for x in row]) + '\n')
        write_row(['tax_name', 'benchmark', 'gu_pairing', 'mismatches',
            'value'])
        for r in benchmark_results:
            write_row(r)


def run(kmer_sample_frac, kmers28, tax_ids,
        tax_ids_to_sample=None):
    # Subsample k-mers
    if kmer_sample_frac is not None:
        num_to_insert = max(1, int(len(kmers28) * kmer_sample_frac))
        kmers_to_insert = set(random.sample(kmers28.keys(), num_to_insert))
        kmers28_subsampled = {}
        for kmer in kmers_to_insert:
            kmers28_subsampled[kmer] = kmers28[kmer]
        kmers28 = kmers28_subsampled

    logging.info("Building trie of 28-mers, forward")
    t_0 = benchmark.build_trie(kmers28, kmer_sample_frac=None, rm=False)

    logging.info("Building trie of 28-mers, reverse")
    # Insert the reverse of each k-mer
    kmers28_rev = {}
    for kmer in kmers28:
        kmer_rev = kmer[::-1]
        kmers28_rev[kmer_rev] = kmers28[kmer]
    t_1 = benchmark.build_trie(kmers28_rev, kmer_sample_frac=None, rm=True)
    del kmers28
    del kmers28_rev

    # Benchmark for each taxonomy
    logging.info("Benchmarking queries on the tries across %d taxonomies",
            len(tax_ids))
    benchmark_queries_across_taxonomies(t_0, t_1, tax_ids,
            ('out/benchmark-queries-pigeonhole-mismatches-subsample-' +
            str(kmer_sample_frac) + '.tsv'),
            tax_ids_to_sample=tax_ids_to_sample)


def main():
    logging.info("Reading sequences")
    seqs = read_kmers.read_seqs('../data/fastas')

    # Build the trie of 28-kmers
    logging.info("Parsing 28-mers from sequences")
    kmers28, tax_ids = read_kmers.read_kmers(seqs, k=28)

    for kmer_sample_frac in [0.0001, 0.0002, 0.0004, 0.0008, 0.0016, 0.0032]:
       if kmer_sample_frac is None:
           # Using all kmers, so only query for 10 of the tax ids (so
           # it's not so slow
           s = 10
       else:
           # Using just a fraction of kmers, so query for all tax ids
           s = None
       run(kmer_sample_frac, kmers28, tax_ids, tax_ids_to_sample=s)


if __name__ == "__main__":
    main()
