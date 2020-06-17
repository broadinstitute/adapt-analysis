"""Benchmark trie data structure, with queries based on partitioning.

The basic idea here is: instead of building a trie of 28-mers and querying
a full 28-mer up to 4 mismatches, build a trie of 14-mers and query
each 14-mer in a 28-mer up to 2 mismatches and return it is present iff
at least one of the two has a hit. This permits false positives but not
false negatives.
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


def query_for_taxonomy(t_14, t_28, taxid, kmer_sample_size=100):
    """Query random sampling of k-mers from a given taxonomy in the trie.

    This first masks that taxonomy from the trie (otherwise most results
    will be for that taxonomy).

    It splits the k-mer (assuming it is a 28-mer) into two 14-mer partitions.
    It then queries each of these

    Args:
        t_14: trie.Trie object containing 14-mers
        t_28: trie.Trie object containing 28-mers
        taxid: taxonomy identifier
        kmer_sample_size: number of k-mers to sample for querying

    Returns:
        (num_matches, num_nodes_visited, runtime, num_matches_from_28_mer) where
        each is a dict {(gu_pairing, mismatches): v} where v is a list of
        values across the randomly sampled k-mers
    """
    # Find all the k-mers in taxid
    logging.info("Finding k-mers for taxid %d", taxid)
    leaves = t_28.root_node.traverse_and_find(taxid)

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
        return None

    # Randomly sample k-mers with replacement, weighting the selection by
    # the number of sequences in taxid that contain each k-mer
    taxid_kmers_sample = random.choices(taxid_kmers,
            weights=taxid_kmer_occ,
            k=kmer_sample_size)

    # Mask taxid from the tries
    logging.info("Masking taxid %d", taxid)
    t_14.mask(taxid)
    t_28.mask(taxid)

    logging.info("Querying for taxid %d", taxid)

    # Query each k-mer with different parameters
    num_matches_kmer1 = {}
    num_matches_kmer2 = {}
    num_matches_combined = {}
    perf_num_nodes_visited = {}
    perf_runtime = {}
    num_matches_from_28_mer = {}
    num_nodes_visited_from_28_mer = {}
    for gu_pairing in [False, True]:
        for m in [0, 1, 2, 3, 4]:
            logging.info("Querying with GU-pairing=%s and mismatches=%d",
                    str(gu_pairing), m)
            num_matches_kmer1_for_params = []
            num_matches_kmer2_for_params = []
            num_matches_combined_for_params = []
            perf_num_nodes_visited_for_params = []
            perf_runtime_for_params = []
            num_matches_from_28_mer_for_params = []
            num_nodes_visited_from_28_mer_for_params = []
            for kmer in taxid_kmers_sample:
                # Query based on a partition into 14-mers, using 1/2
                # (floor'd) the number of mismatches
                start_time = time.time()
                m_half = int(m / 2.0)
                assert len(kmer) == 28
                kmer_1 = kmer[:14]
                kmer_2 = kmer[14:]
                kmer_1_results, kmer_1_num_nodes_visited = t_14.query(
                        kmer_1, mismatches=m_half, gu_pairing=gu_pairing)
                kmer_2_results, kmer_2_num_nodes_visited = t_14.query(
                        kmer_2, mismatches=m_half, gu_pairing=gu_pairing)
                combined_results = benchmark.KmerLeaf([])
                for r in kmer_1_results + kmer_2_results:
                    combined_results.extend(r)
                end_time = time.time()

                num_matches_kmer1_for_params += [len(kmer_1_results)]
                num_matches_kmer2_for_params += [len(kmer_2_results)]
                num_matches_combined_for_params += [len(kmer_1_results) +
                        len(kmer_2_results)]
                num_nodes_visited_combined = kmer_1_num_nodes_visited + kmer_2_num_nodes_visited
                perf_num_nodes_visited_for_params += [num_nodes_visited_combined]
                perf_runtime_for_params += [end_time - start_time]

                # Also check for k-mer in the trie of full 28-mers, to
                # determine if any results are false positives
                full_kmer_results, full_kmer_num_nodes_visited = t_28.query(
                        kmer, mismatches=m, gu_pairing=gu_pairing)
                num_matches_from_28_mer_for_params += [len(full_kmer_results)]
                num_nodes_visited_from_28_mer_for_params += [full_kmer_num_nodes_visited]
            num_matches_kmer1[(gu_pairing, m)] = num_matches_kmer1_for_params
            num_matches_kmer2[(gu_pairing, m)] = num_matches_kmer2_for_params
            num_matches_combined[(gu_pairing, m)] = num_matches_combined_for_params
            perf_num_nodes_visited[(gu_pairing, m)] = perf_num_nodes_visited_for_params
            perf_runtime[(gu_pairing, m)] = perf_runtime_for_params
            num_matches_from_28_mer[(gu_pairing, m)] = num_matches_from_28_mer_for_params
            num_nodes_visited_from_28_mer[(gu_pairing, m)] = num_nodes_visited_from_28_mer_for_params

    # Unmask taxid
    logging.info("Unmasking taxid %d", taxid)
    t_14.unmask_all()
    t_28.unmask_all()

    return (num_matches_kmer1, num_matches_kmer2, num_matches_combined,
            perf_num_nodes_visited, perf_runtime,
            num_matches_from_28_mer, num_nodes_visited_from_28_mer)


def benchmark_queries_across_taxonomies(t_14, t_28, tax_ids, out_tsv):
    """Benchmark queries across taxonomies.

    Args:
        t_14: trie.Trie object containing 14-mers
        t_28: trie.Trie object containing 28-mers
        tax_ids: dict {tax_name: taxonomy identifier} where taxonomy
            identifier is simply an integer
        out_tsv: path to TSV file to which to write benchmark results
    """
    benchmark_results = []
    for tax_name, taxid in tax_ids.items():
        x = query_for_taxonomy(t_14, t_28, taxid)
        if x is None:
            # No k-mers for taxid
            continue
        num_matches_kmer1, num_matches_kmer2, num_matches_combined, perf_num_nodes_visited, perf_runtime, num_matches_from_28_mer, num_nodes_visited_from_28_mer = x
        for benchmark_name, d in [('matches_kmer1', num_matches_kmer1),
                ('matches_kmer2', num_matches_kmer2),
                ('matches_combined', num_matches_combined),
                ('nodes_visited', perf_num_nodes_visited),
                ('runtime', perf_runtime),
                ('matches_from_28_mer', num_matches_from_28_mer),
                ('nodes_visited_from_28_mer', num_nodes_visited_from_28_mer)]:
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


def run(kmer_sample_frac, kmers28, tax_ids):
    # Subsample k-mers before making 14-mers; this is better than
    # passing kmer_sample_frac to benchmark.build_trie() because it
    # ensures 14-mers come from the same subsample of 28-mers
    num_to_insert = max(1, int(len(kmers28) * kmer_sample_frac))
    kmers_to_insert = set(random.sample(kmers28.keys(), num_to_insert))
    kmers28_subsampled = {}
    for kmer in kmers_to_insert:
        if kmer in kmers_to_insert:
            kmers28_subsampled[kmer] = kmers28[kmer]
    kmers28 = kmers28_subsampled

    # Build the trie of 14-kmers
    # Use 14-mers from the same 28-mers selected above
    logging.info("Parsing 14-mers from sequences")
    kmers14 = defaultdict(set)
    for kmer28, tax_with_kmer in kmers28.items():
        for kmer14 in [kmer28[i:(i+14)] for i in range(len(kmer28) - 14 + 1)]:
            for v in tax_with_kmer:
                kmers14[kmer14].add(v)

    logging.info("Building trie of 28-mers")
    t_28 = benchmark.build_trie(kmers28, kmer_sample_frac=None)
    del kmers28

    logging.info("Building trie of 14-mers")
    t_14 = benchmark.build_trie(kmers14, kmer_sample_frac=None)
    del kmers14

    # Benchmark for each taxonomy
    logging.info("Benchmarking queries on the trie across %d taxonomies",
            len(tax_ids))
    benchmark_queries_across_taxonomies(t_14, t_28, tax_ids,
            ('out/benchmark-queries-partition-14mers-subsample-' +
            str(kmer_sample_frac) + '.tsv'))


def main():
    logging.info("Reading sequences")
    seqs = read_kmers.read_seqs('../data/fastas')

    # Build the trie of 28-kmers
    logging.info("Parsing 28-mers from sequences")
    kmers28, tax_ids = read_kmers.read_kmers(seqs, k=28)

    for kmer_sample_frac in [0.0001, 0.0002, 0.0004, 0.0008, 0.0016, 0.0032,
            0.0064, 0.0128]:
        run(kmer_sample_frac, kmers28, tax_ids)


if __name__ == "__main__":
    main()
