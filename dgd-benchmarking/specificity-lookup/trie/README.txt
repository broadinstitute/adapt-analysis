The key results are produced by `benchmark_trie_signatures.py` and plotted by
`plot_benchmark_query_results_for_signature_scheme.R`. This 'signature' scheme
(aka the 'sharding' scheme) is described at the top of
`benchmark_trie_signatures.py`: it produces signatures from 14 or 28-mers and
shards the 28-mers based on the signatures. The same Python script also records
and reports benchmark results on a single, regular trie (no sharding), and the
R script plots these results as well.

Note that `benchmark.py` benchmarks on a single, regular trie, but the
advantage to doing it in `benchmark_trie_signatures.py` as well is that it
allows benchmarking on the same data (taxonomies and k-mers) as the sharding
approaches.

`benchmark_pigeonhole_mismatches.py` is useful for seeing how well one aspect
of the sharding approach works: only querying up to floor(m/2) mismatches for
each half. It appears this does help, but only provides part of the benefit of
the sharding approach.

`benchmark_partition.py` implements a scheme that did not pan out. See emails
for details (in short, it leads to too many false positives; i.e., only a
tiny fraction of the queries it says are non-specific are truly non-specific).
