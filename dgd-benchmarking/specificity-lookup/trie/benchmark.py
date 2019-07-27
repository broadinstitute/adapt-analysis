"""Benchmark trie data structure.
"""

import logging
import random
import time

import read_kmers
import trie

__author__ = 'Hayden Metsky <hayden@mit.edu>'

# Configure basic logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)


class KmerLeaf(trie.LeafInfo):
    """Store info in each leaf node on the k-mer leading to the leaf.

    This stores a dict, d: {taxonomy identifier: {sequence identifiers}}
    giving the taxonomies and sequences containing the k-mer.
    """
    def __init__(self, xs):
        """
        Args:
            xs: collection of tuples (a, b) where a gives taxonomy identifier
                and b gives sequence identifier for a sequence containing
                this k-mer
        """
        self.d = {}
        for x in xs:
            taxid, seqid = x
            if taxid not in self.d:
                self.d[taxid] = set()
            self.d[taxid].add(seqid)

    def extend(self, other):
        """Extend self to contain what is in other.

        Args:
            other: KmerLeaf object
        """
        for taxid, seqids in other.d.items():
            if taxid not in self.d:
                self.d[taxid] = set()
            self.d[taxid].update(seqids)

    def remove(self, taxid):
        """Remove taxonomy identifier from self.

        Args:
            taxid: taxonomy identifier
        """
        if taxid in self.d:
            del self.d[taxid]

    def is_empty(self):
        """Determine if self has any information stored.

        Returns:
            True iff self.d is empty
        """
        return len(self.d) == 0

    def __contains__(self, taxid):
        """Determine if a taxonomy identifier is stored by this leaf.

        Args:
            taxid: taxonomy identifier

        Returns:
            True iff taxid is in self.d
        """
        return taxid in self.d

    def copy(self):
        """Copy self.

        Returns:
            KmerLeaf object, identical to this one
        """
        d_new = {}
        for taxid in self.d.keys():
            d_new[taxid] = set(self.d[taxid])
        new = KmerLeaf({(-1,-1)})
        new.d = d_new
        return new


def build_trie(kmers, kmer_sample_frac=0.0128):
    """Build a trie from 28-mers.

    While building, this simultaenously removes k-mers from the input kmers
    dict so that both the trie and kmers dict do not need to be stored
    at the same time (both may be memory-intensive).

    Args:
        kmers: dict {kmer: {(taxonomy identifier, sequence id)}}
        kmer_sample_frac: fraction of all k-mers to insert into the tree
            (randomly sampled); if None, use all

    Returns:
        trie.Trie object
    """
    if kmer_sample_frac is None:
        # Iterate over set(kmers.keys()) so that dict (kmers) does not change
        # size during iteration
        kmers_to_insert = set(kmers.keys())
    else:
        num_to_insert = max(1, int(len(kmers) * kmer_sample_frac))
        kmers_to_insert = set(random.sample(kmers.keys(), num_to_insert))

    t = trie.Trie()
    i = 1
    num_kmers = len(kmers_to_insert)
    for kmer in kmers_to_insert:
        if i % 10000 == 0:
            logging.info("Inserting k-mer %d of %d", i, num_kmers)
        i += 1

        leaf = KmerLeaf(kmers[kmer])
        t.insert([(kmer, leaf)])
        del kmers[kmer]
    return t


def build_trie_space(kmers, rm=False):
    """Build a space of tries from 28-mers.

    Args:
        kmers_with_sig: dict {kmer: (signature, {(taxonomy identifier, sequence id)})}
        rm: if True, remove k-mers from the input kmers dict while building
            so that they both do not need to be stored at the same time

    Returns:
        trie.TrieSpace object
    """
    # Iterate over set(kmers.keys()) so that dict (kmers) does not change
    # size during iteration
    kmers_to_insert = set(kmers.keys())

    ts = trie.TrieSpace()
    i = 1
    num_kmers = len(kmers_to_insert)
    for kmer in kmers_to_insert:
        if i % 10000 == 0:
            logging.info("Inserting k-mer %d of %d", i, num_kmers)
        i += 1

        sig, vals = kmers[kmer]
        t = ts.get(sig)
        leaf = KmerLeaf(vals)
        t.insert([(kmer, leaf)])
        if rm:
            del kmers[kmer]
    return ts


def query_for_taxonomy(t, taxid, kmer_sample_size=100):
    """Query random sampling of k-mers from a given taxonomy in the trie.

    This first masks that taxonomy from the trie (otherwise most results
    will be for that taxonomy).

    Args:
        t: trie.Trie object
        taxid: taxonomy identifier
        kmer_sample_size: number of k-mers to sample for querying

    Returns:
        (num_matches, num_nodes_visited, runtime) where each is a dict
        {(gu_pairing, mismatches): v} where v is a list of values across
        the randomly sampled k-mers
    """
    # Find all the k-mers in taxid
    logging.info("Finding k-mers for taxid %d", taxid)
    leaves = t.root_node.traverse_and_find(taxid)

    # Construct a list of k-mers and a corresponding list of the number of
    # sequences in taxid that contain each k-mer, to use for weighting
    # during sampling
    taxid_kmers = []
    taxid_kmer_occ = []
    for kmer, li in leaves:
        num_seqs_with_kmer_in_taxid = len(li.d[taxid])
        taxid_kmers += [kmer]
        taxid_kmer_occ += [num_seqs_with_kmer_in_taxid]

    # Randomly sample k-mers with replacement, weighting the selection by
    # the number of sequences in taxid that contain each k-mer
    taxid_kmers_sample = random.choices(taxid_kmers,
            weights=taxid_kmer_occ,
            k=kmer_sample_size)

    # Mask taxid from the trie
    logging.info("Masking taxid %d", taxid)
    t.mask(taxid)

    logging.info("Querying for taxid %d", taxid)

    # Query each k-mer with different parameters
    num_matches = {}
    perf_num_nodes_visited = {}
    perf_runtime = {}
    for gu_pairing in [False, True]:
        for m in [0, 1, 2, 3, 4]:
            logging.info("Querying with GU-pairing=%s and mismatches=%d",
                    str(gu_pairing), m)
            num_matches_for_params = []
            perf_num_nodes_visited_for_params = []
            perf_runtime_for_params = []
            for kmer in taxid_kmers_sample:
                start_time = time.time()
                results, num_nodes_visited = t.query(
                        kmer, mismatches=m, gu_pairing=gu_pairing)
                end_time = time.time()
                num_matches_for_params += [len(results)]
                perf_num_nodes_visited_for_params += [num_nodes_visited]
                perf_runtime_for_params += [end_time - start_time]
            num_matches[(gu_pairing, m)] = num_matches_for_params
            perf_num_nodes_visited[(gu_pairing, m)] = perf_num_nodes_visited_for_params
            perf_runtime[(gu_pairing, m)] = perf_runtime_for_params

    # Unmask taxid
    logging.info("Unmasking taxid %d", taxid)
    t.unmask_all()

    return (num_matches, perf_num_nodes_visited, perf_runtime)


def benchmark_queries_across_taxonomies(t, tax_ids, out_tsv):
    """Benchmark queries across taxonomies.

    Args:
        t: trie.Trie object
        tax_ids: dict {tax_name: taxonomy identifier} where taxonomy
            identifier is simply an integer
        out_tsv: path to TSV file to which to write benchmark results
    """
    benchmark_results = []
    for tax_name, taxid in tax_ids.items():
        num_matches, perf_num_nodes_visited, perf_runtime = query_for_taxonomy(t, taxid)
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


def main():
    logging.info("Reading sequences")
    seqs = read_kmers.read_seqs('../data/fastas')

    logging.info("Computing per-taxonomy stats")
    read_kmers.compute_taxonomy_stats(seqs, 'out/stats.tax.tsv')

    logging.info("Parsing k-mers from sequences")
    kmers, tax_ids = read_kmers.read_kmers(seqs, k=28)

    # Compute distribution of number of taxonomies each k-mer appears in
    logging.info("Computing distribution of k-mer occurrences")
    read_kmers.compute_kmer_occurrence_distribution(kmers, 'out/stats.kmer-occ.tsv')

    # Build the trie
    logging.info("Building trie")
    kmer_sample_frac = 0.0128
    t = build_trie(kmers, kmer_sample_frac=kmer_sample_frac)
    del kmers
    logging.info("Done building trie")

    # Benchmark for each taxonomy
    logging.info("Benchmarking queries on the trie across %d taxonomies",
            len(tax_ids))
    benchmark_queries_across_taxonomies(t, tax_ids,
            'out/benchmark-queries-subsample-' + str(kmer_sample_frac) + '.tsv')


if __name__ == "__main__":
    main()
