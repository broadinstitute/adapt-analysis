This analysis measures variability in design from the same underlying input
set of data: with randomly resampled input and also by re-running on the same
full set of data.

Run:
  1) `./make_designs.sh` to produce all designs that will be used for evaluating dispersion
        - This designs in two ways: using randomly resampled input (bootstrapping) and
          using all accessions (to determine the dispersion due to randomness in the design,
          but with the same exact input; non-resampled)

Measuring dispersion:
  2) `./evaluate_dispersion.sh` to use those designs to output a summary of the dispersion across designs
        - This evaluates dispersion for designs both on the resampled and non-resampled inputs

Measuring coverage of resampled designs against all data:
  3) `./compute_coverage_on_all_data.sh` to compute the coverage against all input data
        - This takes only the designs on the resampled inputs, and computes the coverage they
          achieve against all sequences (the non-resampled designs should, inherent in the design,
          achieve the desired coverage against all sequences because they were designed using
          all sequences)
