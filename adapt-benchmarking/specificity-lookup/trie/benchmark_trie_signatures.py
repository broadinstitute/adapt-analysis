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
and the same for q_2. The results with this are the 'split' approach
('split_approach_*').

Each trie should have, on average, 1/2^14 of the total number of 28-mers. So
the total space of 28-mers that we're querying is 30/2^14 = 0.002.


This also benchmarks a related, simplified approach that is based on just
one signature for each 28-mer. There are 2^28 tries each with, on average,
1/2^28 of the total number of 28-mers. To query up to 3 mismatches, we have
to flip the letter at every combination of 1, 2, and 3 positions of the
signature. The results on this are the 'full' approach ('full_approach_*').

Optionally, this will benchmark results on a regular/simple trie (i.e.,
without the sharding approach described above). The results on this are
the 'noshard' approach ('noshard_approach_*').
"""

from collections import defaultdict
from itertools import combinations
import logging
import math
import random
import time

import benchmark
import read_kmers
import trie

__author__ = 'Hayden Metsky <hayden@mit.edu>'

# Set whether to verify all results using a regular/simple trie, and to
# report benchmarking results on this
VERIFY = True

# Configure basic logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


def signatures_with_mismatches(sig, m):
    """Generate signatures within m mismatches given a signature.

    Args:
        sig: signature as string
        m: number of mismatches such that all signatures have 0,1,...,m
            mismatches relative to kmer

    Returns:
        list of signatures as strings
    """
    sigs_with_mismatches = [sig]
    for mi in range(1, m+1):
        # Introduce mi mismatches to sig
        for pos in combinations(range(len(sig)), mi):
            # Introduce mismatches to sig at positions in pos
            # Note that len(pos) == mi
            sig_mi = list(str(sig))
            for j in pos:
                # Introduce a mismatch at position j
                if sig_mi[j] == 'G':
                    sig_mi[j] = 'T'
                else:
                    sig_mi[j] = 'G'
            sigs_with_mismatches += [''.join(sig_mi)]
    return sigs_with_mismatches


def split_signature(kmer, half):
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


def split_signatures_with_mismatches(kmer, half, m):
    """Generate signatures within one mismatch given a k-mer.

    Args:
        kmer: 28-mer
        half: 0 or 1, denoting the first or second half
        m: number of mismatches such that all signatures have 0,1,...,m
            mismatches relative to kmer

    Returns:
        list of signatures as strings
    """
    sig = split_signature(kmer, half)
    return signatures_with_mismatches(sig, m)


def full_signature(kmer):
    """Generate a 2-letter alphabet signature of a full k-mer.

    The transformations are A->G and C->T.

    Args:
        kmer: 28-mer

    Returns:
        signature as a string
    """
    assert len(kmer) == 28
    return kmer.replace('A', 'G').replace('C', 'T')


def full_signatures_with_mismatches(kmer, m):
    """Generate signatures within m mismatches given a k-mer.

    Args:
        kmer: 28-mer
        m: number of mismatches such that all signatures have 0,1,...,m
            mismatches relative to kmer

    Returns:
        list of signatures as strings
    """
    sig = full_signature(kmer)
    return signatures_with_mismatches(sig, m)


def query_for_taxonomy(ts_0, ts_1, ts_full, taxid, kmer_sample_size=100,
        t_for_verification=None):
    """Query random sampling of k-mers from a given taxonomy in the trie
    spaces.

    This first masks that taxonomy from the tries (otherwise most results
    will be for that taxonomy).

    Args:
        ts_0: trie.TrieSpace object key'd on first 14-mer
        ts_1: trie.TrieSpace object key'd on first 14-mer
        ts_full: trie.TrieSpace object key'd on full 28-mer
        taxid: taxonomy identifier
        kmer_sample_size: number of k-mers to sample for querying
        t_for_verification: if set, a basic trie with all 28-mers to use
            for verifying query results and reporting benchmarks on this
            ('noshard' results)

    Returns:
        tuple (split, full) where each is a tuple giving information on the
        split signature approach or full signature approach:
            (has_hit, num_nodes_visited, runtime) where
            each is a dict {(gu_pairing, mismatches): v} where v is a list of
            values across the randomly sampled k-mers
    """
    # Find all the k-mers in taxid
    # Only do this using the split trie spaces (ts_0 and ts_1); ts_full should
    # provide the same k-mers
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
        logging.info("No k-mers for taxid %d", taxid)
        return None

    # Randomly sample k-mers with replacement, weighting the selection by
    # the number of sequences in taxid that contain each k-mer
    taxid_kmers_sample = random.choices(taxid_kmers,
            weights=taxid_kmer_occ,
            k=kmer_sample_size)

    # Mask taxid from the tries
    logging.info("Masking taxid %d", taxid)
    for ts in [ts_0, ts_1, ts_full]:
        for sig, t in ts.all_tries():
            t.mask(taxid)
    if t_for_verification is not None:
        t_for_verification.mask(taxid)

    logging.info("Querying for taxid %d", taxid)

    # Query each k-mer at different parameters (but all with G-U pairing)
    split_has_hit = {}
    split_num_results = {}
    split_perf_num_nodes_visited = {}
    split_perf_runtime = {}
    full_has_hit = {}
    full_num_results = {}
    full_perf_num_nodes_visited = {}
    full_perf_runtime = {}
    noshard_has_hit = {}
    noshard_num_results = {}
    noshard_perf_num_nodes_visited = {}
    noshard_perf_runtime = {}
    def evaluate_for_params(gu_pairing, m):
        split_has_hit_for_params = []
        split_num_results_for_params = []
        split_perf_num_nodes_visited_for_params = []
        split_perf_runtime_for_params = []
        full_has_hit_for_params = []
        full_num_results_for_params = []
        full_perf_num_nodes_visited_for_params = []
        full_perf_runtime_for_params = []
        noshard_has_hit_for_params = []
        noshard_num_results_for_params = []
        noshard_perf_num_nodes_visited_for_params = []
        noshard_perf_runtime_for_params = []
        for kmer in taxid_kmers_sample:
            ##########
            # Benchmark the approach that splits a 28-mer into two halves
            # and computes a signature on each half

            start_time = time.time()

            # By the pigeonhole principle, we'll have at most floor(m/2)
            # mismatches in one of the two halves of the 28-mer
            m_half = int(m / 2)

            # Determine signatures of each half of the k-mer to use
            if m_half == 0:
                sigs1 = [split_signature(kmer, 0)]
                sigs2 = [split_signature(kmer, 1)]
            else:
                sigs1 = split_signatures_with_mismatches(kmer, 0, m_half)
                sigs2 = split_signatures_with_mismatches(kmer, 1, m_half)

            num_nodes_visited = 0
            all_results = []

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
                        num_nodes_visited += q_num_nodes_visited
                        all_results.extend(q_results)

            end_time = time.time()

            # Just adding up len(q_results) above can overcount because the
            # same result (hit) can appear for queries from both the first and
            # second half of the kmer; instead, simply compute the length
            # of the union of all_results
            results_of_queries = benchmark.KmerLeaf.union(all_results)
            num_results = len(results_of_queries)
            has_hit = int(num_results > 0)

            split_has_hit_for_params += [has_hit]
            split_num_results_for_params += [num_results]
            split_perf_num_nodes_visited_for_params += [num_nodes_visited]
            split_perf_runtime_for_params += [end_time - start_time]
            ##########

            ##########
            # Verify the results and benchmark the simple/no-sharding approach
            if t_for_verification is not None:
                start_time = time.time()
                true_results, num_nodes_visited = t_for_verification.query(
                        kmer, mismatches=m, gu_pairing=gu_pairing)
                end_time = time.time()

                true_results = benchmark.KmerLeaf.union(true_results)
                if results_of_queries != true_results:
                    raise Exception(("Trie verification failed (split approach)"))

                num_results = len(true_results)
                has_hit = int(num_results > 0)

                noshard_has_hit_for_params += [has_hit]
                noshard_num_results_for_params += [num_results]
                noshard_perf_num_nodes_visited_for_params += [num_nodes_visited]
                noshard_perf_runtime_for_params += [end_time - start_time]
            ##########

            ##########
            # Benchmark the related approach of using a signature based on
            # a full 28-mer

            start_time = time.time()

            all_results = []
            num_nodes_visited = 0

            sigs = full_signatures_with_mismatches(kmer, m)
            for sig in sigs:
                t = ts_full.get(sig, make=False)
                if t is None:
                    # No results
                    pass
                else:
                    q = kmer
                    q_results, q_num_nodes_visited = t.query(
                            q, mismatches=m, gu_pairing=gu_pairing)
                    num_nodes_visited += q_num_nodes_visited
                    all_results.extend(q_results)

            end_time = time.time()

            results_of_queries = benchmark.KmerLeaf.union(all_results)
            num_results = len(results_of_queries)
            has_hit = int(num_results > 0)

            # Verify the results
            if t_for_verification is not None:
                if results_of_queries != true_results:
                    raise Exception(("Trie verification failed (full approach)"))

            full_has_hit_for_params += [has_hit]
            full_num_results_for_params += [num_results]
            full_perf_num_nodes_visited_for_params += [num_nodes_visited]
            full_perf_runtime_for_params += [end_time - start_time]
            ##########

        split_has_hit[(gu_pairing, m)] = split_has_hit_for_params
        split_num_results[(gu_pairing, m)] = split_num_results_for_params
        split_perf_num_nodes_visited[(gu_pairing, m)] = split_perf_num_nodes_visited_for_params
        split_perf_runtime[(gu_pairing, m)] = split_perf_runtime_for_params
        full_has_hit[(gu_pairing, m)] = full_has_hit_for_params
        full_num_results[(gu_pairing, m)] = full_num_results_for_params
        full_perf_num_nodes_visited[(gu_pairing, m)] = full_perf_num_nodes_visited_for_params
        full_perf_runtime[(gu_pairing, m)] = full_perf_runtime_for_params
        noshard_has_hit[(gu_pairing, m)] = noshard_has_hit_for_params
        noshard_num_results[(gu_pairing, m)] = noshard_num_results_for_params
        noshard_perf_num_nodes_visited[(gu_pairing, m)] = noshard_perf_num_nodes_visited_for_params
        noshard_perf_runtime[(gu_pairing, m)] = noshard_perf_runtime_for_params

    for gu_pairing in [False, True]:
        for m in [0, 1, 2, 3, 4, 5, 6]:
            evaluate_for_params(gu_pairing, m)

    # Unmask taxid
    logging.info("Unmasking taxid %d", taxid)
    for ts in [ts_0, ts_1, ts_full]:
        for sig, t in ts.all_tries():
            t.unmask_all()
    if t_for_verification is not None:
        t_for_verification.unmask_all()

    return ((split_has_hit, split_num_results,
                split_perf_num_nodes_visited, split_perf_runtime),
            (full_has_hit, full_num_results,
                full_perf_num_nodes_visited, full_perf_runtime),
            (noshard_has_hit, noshard_num_results,
                noshard_perf_num_nodes_visited, noshard_perf_runtime))


def benchmark_queries_across_taxonomies(ts_0, ts_1, ts_full, tax_ids, out_tsv,
        tax_ids_to_sample=None, t_for_verification=None):
    """Benchmark queries across taxonomies.

    Args:
        ts_0: trie.TrieSpace object key'd on first 14-mer
        ts_1: trie.TrieSpace object key'd on first 14-mer
        ts_full: trie.TrieSpace object key'd on full 28-mer
        tax_ids: dict {tax_name: taxonomy identifier} where taxonomy
            identifier is simply an integer
        out_tsv: path to TSV file to which to write benchmark results
        tax_ids_to_sample: number of tax ids to sample
        t_for_verification: if set, a basic trie with all 28-mers to use
            for verifying query results and benchmarking
    """
    if tax_ids_to_sample is not None:
        tax_ids_to_insert = set(random.sample(tax_ids.keys(), tax_ids_to_sample))
        tax_ids_subsampled = {}
        for k in tax_ids_to_insert:
            tax_ids_subsampled[k] = tax_ids[k]
        tax_ids = tax_ids_subsampled

    benchmark_results = []
    for tax_name, taxid in tax_ids.items():
        x = query_for_taxonomy(ts_0, ts_1, ts_full, taxid,
                t_for_verification=t_for_verification)
        if x is None:
            # No k-mers for taxid
            continue
        split_info, full_info, noshard_info = x
        split_has_hit, split_num_results, split_perf_num_nodes_visited, split_perf_runtime = split_info
        full_has_hit, full_num_results, full_perf_num_nodes_visited, full_perf_runtime = full_info
        noshard_has_hit, noshard_num_results, noshard_perf_num_nodes_visited, noshard_perf_runtime = noshard_info
        for benchmark_name, d in [('split_approach_has_hit', split_has_hit),
                ('split_approach_num_results', split_num_results),
                ('split_approach_nodes_visited', split_perf_num_nodes_visited),
                ('split_approach_runtime', split_perf_runtime),
                ('full_approach_has_hit', full_has_hit),
                ('full_approach_num_results', full_num_results),
                ('full_approach_nodes_visited', full_perf_num_nodes_visited),
                ('full_approach_runtime', full_perf_runtime),
                ('noshard_approach_has_hit', noshard_has_hit),
                ('noshard_approach_num_results', noshard_num_results),
                ('noshard_approach_nodes_visited', noshard_perf_num_nodes_visited),
                ('noshard_approach_runtime', noshard_perf_runtime)]:
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
        tax_ids_to_sample=None, verify_in_trie=False):
    # Subsample k-mers
    if kmer_sample_frac is not None:
        num_to_insert = max(1, int(len(kmers28) * kmer_sample_frac))
        kmers_to_insert = set(random.sample(kmers28.keys(), num_to_insert))
        kmers28_subsampled = {}
        for kmer in kmers_to_insert:
            kmers28_subsampled[kmer] = kmers28[kmer]
        kmers28 = kmers28_subsampled

    if verify_in_trie:
        # Pass the same already-downsampled kmers (kmers28) rather than
        # downsampling separately
        logging.info("Building basic trie for verifying results")
        t_for_verification = benchmark.build_trie(
                kmers28,
                kmer_sample_frac=None,
                rm=False)
    else:
        t_for_verification = None

    logging.info("Building trie space of 28-mers, key'd on the full 28-mer")
    # Build trie space key'd by the full 28-mer, for the approach that
    # uses the full signature (rather than splitting into two halves)
    for kmer in kmers28:
        sig = full_signature(kmer)
        kmers28[kmer] = (sig, kmers28[kmer])
    ts_full = benchmark.build_trie_space(kmers28, rm=False)

    logging.info("Building trie space of 28-mers, key'd by the first 14-mer")
    # Build trie space key'd by the first 14-mer
    for kmer in kmers28:
        sig = split_signature(kmer, 0)
        old_sig, vals = kmers28[kmer]
        kmers28[kmer] = (sig, vals)
    ts_0 = benchmark.build_trie_space(kmers28, rm=False)

    logging.info("Building trie space of 28-mers, key'd by the second 14-mer")
    # Build trie space key'd by the second 14-mer; put kmers into this trie
    # in reverse, to make queries faster
    kmers28_rev = {}
    for kmer in kmers28:
        sig = split_signature(kmer, 1)
        old_sig, vals = kmers28[kmer]
        kmer_rev = kmer[::-1]
        kmers28_rev[kmer_rev] = (sig, vals)
    ts_1 = benchmark.build_trie_space(kmers28_rev, rm=True)
    del kmers28
    del kmers28_rev

    # Benchmark for each taxonomy
    logging.info("Benchmarking queries on the data structure across %d taxonomies",
            len(tax_ids))
    benchmark_queries_across_taxonomies(ts_0, ts_1, ts_full, tax_ids,
            ('out/benchmark-queries-signature-subsample-' +
            str(kmer_sample_frac) + '.tsv'),
            tax_ids_to_sample=tax_ids_to_sample,
            t_for_verification=t_for_verification)


def main():
    logging.info("Reading sequences")
    seqs = read_kmers.read_seqs('../data/fastas')

    # Build the trie of 28-kmers
    logging.info("Parsing 28-mers from sequences")
    kmers28, tax_ids = read_kmers.read_kmers(seqs, k=28)

    #test_fracs = [0.0001, 0.0002, 0.0004, 0.0008, 0.0016, 0.0032, 0.0128]
    #test_fracs = [0.002]
    #test_fracs = [None]
    test_fracs = [0.0001, 0.0016, 0.0128]
    for kmer_sample_frac in test_fracs:
       if kmer_sample_frac is None:
           # Using all kmers, so only query for 10 of the tax ids (so
           # it's not so slow)
           s = min(10, len(tax_ids))
       elif kmer_sample_frac > 0.01:
           # Use only some tax ids so it's not so slow
           s = min(100, len(tax_ids))
       else:
           # Using just a fraction of kmers, so query for all tax ids
           s = None
       run(kmer_sample_frac, kmers28, tax_ids, tax_ids_to_sample=s,
           verify_in_trie=VERIFY)


if __name__ == "__main__":
    main()
