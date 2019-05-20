"""Benchmark trie data structure.
"""

import logging
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
        for taxid, seqid in other.d.items():
            if taxid not in self.d:
                self.d[taxid] = set()
            self.d[taxid].add(seqid)

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


def build_trie(kmers):
    """Build a trie from 28-mers.

    While building, this simultaenously removes k-mers from the input kmers
    dict so that both the trie and kmers dict do not need to be stored
    at the same time (both may be memory-intensive).

    Args:
        kmers: dict {kmer: {(taxonomy identifier, sequence id)}}

    Returns:
        trie.Trie object
    """
    t = trie.Trie()
    # Iterate over set(kmers.keys()) so that dict (kmers) does not change
    # size during iteration
    i = 1
    num_kmers = len(kmers)
    for kmer in set(kmers.keys()):
        if i % 10000 == 0:
            logging.info("Inserting k-mer %d of %d", i, num_kmers)
        i += 1

        leaf = KmerLeaf(kmers[kmer])
        t.insert([(kmer, leaf)])
        del kmers[kmer]
    return t


def query_for_taxonomy(t, taxid):
    """Query all k-mers from a given taxonomy in the trie.

    This first masks that taxonomy from the trie (otherwise most results
    will be for that taxonomy).

    Args:
        t: trie.Trie object
        taxid: taxonomy identifier
    """
    # Find all the k-mers in taxid
    logging.info("Finding k-mers for taxid %d", taxid)
    leaves = t.root_node.traverse_and_find(taxid)
    taxid_kmers = [kmer for kmer, _ in leaves]

    # Mask taxid from the trie
    logging.info("Masking taxid %d", taxid)
    t.mask(taxid)

    logging.info("Querying for taxid %d", taxid)

    # Query each k-mer
    for gu_pairing in [False, True]:
        for m in [0, 1, 2, 3, 4, 5]:
            total_results = 0
            total_nodes_visited = 0
            start = time.time()
            for kmer in taxid_kmers:
                results, num_nodes_visited = t.query(
                        kmer, mismatches=m, gu_pairing=gu_pairing)
                total_results += len(results)
                total_nodes_visited += num_nodes_visited
            print(gu_pairing, m, total_results, (time.time()-start)/len(taxid_kmers),
                    float(total_nodes_visited)/len(taxid_kmers))

    # Unmask taxid
    logging.info("Unmasking taxid %d", taxid)
    t.unmask_all()


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
    t = build_trie(kmers)
    del kmers
    logging.info("Done building trie")

    for tax_name, taxid in tax_ids.items():
        # Query for taxonomy taxid
        query_for_taxonomy(t, taxid)


if __name__ == "__main__":
    main()
