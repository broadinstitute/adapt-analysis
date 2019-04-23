This analysis measures how the design methods perform against naive design
approaches.

All designs are produced using a sliding window approach, for ease of
comparison. The naive approaches select just one guide per window (e.g.,
the consensus guide that achieves the highest coverage in the window) and
reports on the coverage obtained by the guide in the window (i.e., fraction
of genomes that it covers). The usual design method produces multiple
guides in each window in order to achieve a desired coverage (e.g., 99%
of genomes). As a result, the comparison is not direct: we are comparing
the coverage obtained by a (single guide) naive design against the number
of guides required to achieve a certain coverage.

This analysis is run on multiple inputs, each a random resampling of all
the available genomes. This allows us to construct distributions of values
(e.g., of the coverage obtained by a naive design in each window, or of
the number of guides required to achieve a desired coverage in a window).

Run:
  1) `./make_designs.sh` to produce designs, with both the usual methods and
     naive methods
        - This produces files tax-[taxonomy]/input-alns/design-{i}.fasta.0
          containing an alignment on a resampling of input sequences.
        - This produces files tax-[taxonomy]/designs/design-{i}.real-design.tsv.0,
          where each design-{i}.real-design.tsv.0 gives a design using
          `design.py` with the `sliding-window` approach.
        - This produces files tax-[taxonomy]/designs/design-{i}.naive-design.tsv,
          where each design-{i}.naive-design.tsv gives a design using
          `design_naively.py`

Compile the above results:
  2) `./compile_designs.sh` to compile the outputs of `./make_designs.sh`
        - This produces files tax-[taxonomy]/naive-designs.tsv and
          tax-[taxonomy]/real-designs.tsv, which simply concatenate the outputs
          above across the different {i}
