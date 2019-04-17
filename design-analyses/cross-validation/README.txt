This analysis measures the coverage that design on one set of data achieves on
another set of data.

It is analogous to cross-validation in machine learning. Of all input
accessions, 80% are randomly chosen as a 'design-set' and the other 20% are
selected to be the 'test-set'. Designs are created on the design set, and then
coverage is computed on the test set. This is performed many times to get a
distribution of coverage on the test set.

Run:
  1) `./make_designs_and_test_coverage.sh`
        - This splits the input data into a design set and test set, performs
          the design on the design set, and measures its coverage on the
          test set. It does this many times.

Make a distribution of coverage values:
  2) `./compile_coverages.sh`
        - This reads the coverage values above -- one per target in each
          design -- and summarizes them into a distribution of coverage values.

Plot a distribution of coverage values:
  3) `./make_plots.sh`
        - This plots the distribution of coverage values summarized above.
