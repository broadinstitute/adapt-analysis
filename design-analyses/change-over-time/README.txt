This analysis measures the coverage that designs produced in some year achieve against
sequences collected in a particular year.

Designs 'produced in some year Y' are designed using all available sequences
collected up to and including in year Y. Coverage achieved 'against sequences
collected in a particular year' is measured by seeing what fraction of
sequences collected in that year (and only in that year) are covered by a
design.

This is intended to show, if we were to produce a design for a taxon in year Y,
how that design would hold up over time. It should achieve reasonable coverage
for all years <= Y. (Note, though, that it is possible in theory for it to produce
poor coverage for some of the years <= Y -- e.g., if some year has few sequences
very divergent from all the other years and the desired coverage was <100%. But
the average of sequences covered, weighted by the number of sequences in each year,
should be at least the desired fraction to cover (e.g., >= 95%).) For sequences
collected in years >Y, the coverage achieved may be less than the desired coverage
for the design (and may decline over time).

Run:
  1) `./make_designs.sh` to produce designs using available data up to each year
        - This produces files tax-[taxonomy]/designs/designs_up-to-{Y}/design-{i}.tsv.0,
          where each design-{i}.tsv.0 gives a design using all available sequences
          up to year Y (the different i represent replicates across randomly resampled
          input sequences).
  2) `./compute_coverage_per_year.sh` to determine coverage that the above designs achieve per year
        - This produces files tax-[taxonomy]/designs/designs_up-to-{Y}/coverages/coverage-in-[Z]/design-{i}.coverage.txt,
          which gives the coverage (for each target in the design) that is achieved
          against all sequences in year Z by the design (made with make_designs.sh)
          that was made using all sequences up to year Y. (Again, the different i
          represent replicates.)

 Compile the above results:
  3) `./compile_coverages.sh` compiles the outputs of `./compute_coverage_per_year.sh`
        - This compiles, for each taxonomy, the coverage files produced above into a
          single file. It creates the file tax-[taxonomy]/coverages.distribution.txt.
