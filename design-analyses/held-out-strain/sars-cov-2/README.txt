This evaluates how well designs work on SARS-CoV-2 if made them (on SARS-related CoV species) before we knew about SARS-CoV-2.

See input/README.txt for more details on the input. Note, importantly, that although the
criteria for design input sequences went up to 2019-11-30 (as noted in the filename), the
most recent genome in those sequences is Dec. 31, 2018 -- so it is the same as if our
cutoff were Dec. 31, 2018 or Jan. 1, 2019.

To produce designs:
  1) ./design.sh -- this uses all 311 genomes in the SARS-related CoV species up to 2019-11-30;
     most of them are SARS-CoV-1 (call this design 'A')
     - Output is in designs/design.tsv.0
  2) ./design_downsampled-sars-cov-1.sh -- this is the same as the above, except uses as
     input just a single SARS-CoV-1 genome (49 genomes total) (call this design 'B')
     - Output is in designs/design.downsampled-sars-cov-1.tsv.0
Note that I also made design_minimize-guides.sh, but do not use its output.

To evaluate detection on SARS-CoV-2:
  3) ./evaluate_detection_on_sars-cov-2.sh -- evaluates designs A and B against all
     SARS-CoV-2 genomes available up to 2020-11-12
     - Output for design A is evaluations/design.sars-cov-2.*
     - Output for design B is evaluations/design.downsampled-sars-cov-1.sars-cov-2.*

To evaluate detection against what was used to make the design:
  4) ./evaluate_detection_on_design_input.sh -- evaluates designs A and B against
     the genomes input to the design
     - Output for design A is evaluations/design.design-input.*
     - Output for design B is evaluations/design.downsampled-sars-cov-1.design-input.*

To compile the evaluations into a single file:
  5) ./compile_evaluations.sh -- generates evaluations/evaluations-compiled.tsv, which
     combines the evaluations
