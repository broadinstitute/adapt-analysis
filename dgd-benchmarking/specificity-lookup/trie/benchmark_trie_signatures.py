"""Benchmark trie data structure, with queries based on splitting 28-mers
into separate tries using a 'signature'.

The basic idea here is: for a 28-mer S determine two signatures: s_1 and s_2.
s_1 is the first 14-mer of S collapsed down to a 2-letter alphabet (to
account for G-U pairing). s_2 is the second 14-mer of S collapsed down to a
2-letter alphabet (to account for G-U pairing). We have a separate trie for
every possible signature, and S gets stored twice: once for s_1 and again for
s_2. Now consider a 28-mer query Q, which similarly has signatures q_1 and
q_2. Let's only query for Q up to 3 mismatches. By the pigeonhole principle,
either Q[0:14] or Q[14:28] has <= 1 mismatch with any possible result. So
we can query in 30 tries (corresponding to 30 signatures): q_1, as well as
q_1 where at every position we flip the letter (i.e., introduce a mismatch);
and the same for q_2.

Each trie should have, on average, 1/2^14 of the total number of 28-mers. So
the total space of 28-mers that we're querying is 30/2^14 = 0.002.
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


def signature(kmer, half):
    """Generate a 2-letter alphabet signature of one half of a k-mer.

    The transformations are A->G and C->T.

    Args:
        kmer: 28-mer
        half: 0 or 1, denoting the first or second half

    Returns:
        signature as a string
    """
    assert len(kmer) == 28
    assert half in [0, 1]
    kmer_half = kmer[0:14] if half == 0 else kmer[14:28]
    return kmer_half.replace('A', 'G').replace('C', 'T')


def signatures_with_mismatches(kmer, half):
    """Generate signatures within one mismatch.

    Args:
        kmer: 28-mer
        half: 0 or 1, denoting the first or second half

    Returns:
        list of signatures as strings
    """
    sig = signature(kmer, half)
    sig_one_mismatch = [sig]
    for i in range(len(sig)):
        sig_r = list(str(sig))
        if sig_r[i] == 'G':
            sig_r[i] = 'T'
        else:
            sig_r[i] = 'G'
        sig_one_mismatch += [''.join(sig_r)]
    return sig_one_mismatch


def query_for_taxonomy(ts_0, ts_1, taxid, kmer_sample_size=100):
    """Query random sampling of k-mers from a given taxonomy in the trie
    spaces.

    This first masks that taxonomy from the tries (otherwise most results
    will be for that taxonomy).

    Args:
        t_0: trie.TrieSpace object key'd on first 14-mer
        t_1: trie.TrieSpace object key'd on first 14-mer
        taxid: taxonomy identifier
        kmer_sample_size: number of k-mers to sample for querying

    Returns:
        (num_matches, num_nodes_visited, runtime, num_matches_from_28_mer) where
        each is a dict {(gu_pairing, mismatches): v} where v is a list of
        values across the randomly sampled k-mers
    """
    # Find all the k-mers in taxid
    logging.info("Finding k-mers for taxid %d", taxid)
    leaves = []
    for ts in [ts_0, ts_1]:
        for sig, t in ts.all_tries():
            leaves += t.root_node.traverse_and_find(taxid)

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
    for ts in [ts_0, ts_1]:
        for sig, t in ts.all_tries():
            t.mask(taxid)

    logging.info("Querying for taxid %d", taxid)

    # Query each k-mer at different parameters (but all with G-U pairing)
    num_matches_combined = {}
    perf_num_nodes_visited = {}
    perf_runtime = {}
    gu_pairing = True
    for m in [0, 1, 2, 3]:
        num_matches_combined_for_params = []
        perf_num_nodes_visited_for_params = []
        perf_runtime_for_params = []
        for kmer in taxid_kmers_sample:
            start_time = time.time()

            # By the pigeonhole principle, we'll have at most floor(m/2)
            # mismatches in one of the two halves of the 28-mer
            m_half = int(m / 2)
            assert m_half in [0, 1]

            # Determine signatures of each half of the k-mer to use
            if m_half == 0:
                sigs1 = [signature(kmer, 0)]
                sigs2 = [signature(kmer, 1)]
            else:
                sigs1 = signatures_with_mismatches(kmer, 0)
                sigs2 = signatures_with_mismatches(kmer, 1)

            num_results_combined = 0
            num_nodes_visited_combined = 0

            # Query for tries maybe containing each half of the kmer
            for reverse, sigs, ts in [(False, sigs1, ts_0), (True, sigs2, ts_1)]:
                for sig in sigs:
                    t = ts.get(sig, make=False)
                    if t is None:
                        # No results
                        pass
                    else:
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
                        num_results_combined += len(q_results)
                        num_nodes_visited_combined += q_num_nodes_visited

            end_time = time.time()

            num_matches_combined_for_params += [num_results_combined]
            perf_num_nodes_visited_for_params += [num_nodes_visited_combined]
            perf_runtime_for_params += [end_time - start_time]

        num_matches_combined[(gu_pairing, m)] = num_matches_combined_for_params
        perf_num_nodes_visited[(gu_pairing, m)] = perf_num_nodes_visited_for_params
        perf_runtime[(gu_pairing, m)] = perf_runtime_for_params

    # Unmask taxid
    logging.info("Unmasking taxid %d", taxid)
    for ts in [ts_0, ts_1]:
        for sig, t in ts.all_tries():
            t.unmask_all()

    return (num_matches_combined,
            perf_num_nodes_visited, perf_runtime)


def benchmark_queries_across_taxonomies(t_0, t_1, tax_ids, out_tsv):
    """Benchmark queries across taxonomies.

    Args:
        t_0: trie.TrieSpace object key'd on first 14-mer
        t_1: trie.TrieSpace object key'd on first 14-mer
        tax_ids: dict {tax_name: taxonomy identifier} where taxonomy
            identifier is simply an integer
        out_tsv: path to TSV file to which to write benchmark results
    """
    benchmark_results = []
    for tax_name, taxid in tax_ids.items():
        x = query_for_taxonomy(t_0, t_1, taxid)
        if x is None:
            # No k-mers for taxid
            continue
        num_matches_combined, perf_num_nodes_visited, perf_runtime = x
        for benchmark_name, d in [('matches_combined', num_matches_combined),
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


def run(kmer_sample_frac, kmers28, tax_ids):
    # Subsample k-mers
    num_to_insert = max(1, int(len(kmers28) * kmer_sample_frac))
    kmers_to_insert = set(random.sample(kmers28.keys(), num_to_insert))
    kmers28_subsampled = {}
    for kmer in kmers_to_insert:
        kmers28_subsampled[kmer] = kmers28[kmer]
    kmers28 = kmers28_subsampled

    logging.info("Building trie space of 28-mers, key'd by the first 14-mer")
    # Build trie space key'd by the first 14-mer
    for kmer in kmers28:
        sig = signature(kmer, 0)
        kmers28[kmer] = (sig, kmers28[kmer])
    ts_0 = benchmark.build_trie_space(kmers28, rm=False)

    logging.info("Building trie space of 28-mers, key'd by the second 14-mer")
    # Build trie space key'd by the second 14-mer; put kmers into this trie
    # in reverse, to make queries faster
    kmers28_rev = {}
    for kmer in kmers28:
        sig = signature(kmer, 1)
        old_sig, vals = kmers28[kmer]
        kmer_rev = kmer[::-1]
        kmers28_rev[kmer_rev] = (sig, vals)
    ts_1 = benchmark.build_trie_space(kmers28_rev, rm=True)
    del kmers28
    del kmers28_rev

    # Benchmark for each taxonomy
    logging.info("Benchmarking queries on the trie across %d taxonomies",
            len(tax_ids))
    benchmark_queries_across_taxonomies(ts_0, ts_1, tax_ids,
            ('out/benchmark-queries-signature-subsample-' +
            str(kmer_sample_frac) + '.tsv'))


def main():
    logging.info("Reading sequences")
    seqs = read_kmers.read_seqs('../data/fastas')

    # Build the trie of 28-kmers
    logging.info("Parsing 28-mers from sequences")
    kmers28, tax_ids = read_kmers.read_kmers(seqs, k=28)

    #for kmer_sample_frac in [0.0001, 0.0002, 0.0004, 0.0008, 0.0016, 0.0032,
    #        0.0064, 0.0128]:
    for kmer_sample_frac in [0.0128]:
        run(kmer_sample_frac, kmers28, tax_ids)


if __name__ == "__main__":
    main()
