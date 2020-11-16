This compares runtime with and without having GuideSearcher memoize guide
computations.

The runs without memoization, for some viruses, take too long to complete.
Therefore, rather than simply comparing the total runtime without and
with memoization, we plot the cumulative runtime (only for the search
process - i.e., not end-to-end) vs. the number of windows searched. Then,
at the final window, comparing the runtimes can give a lower bound (assuming
the runtime without memoization *grows* faster than with memoization)
on what the difference in runtime would be.

Steps:
 ( 1) ./run.sh to run ADAPT; the runs without memoization will take too long,
      so this will have to be manually killed at some point
  (2) ./summarize.sh to calculate elapsed runtime at each window during
      the search, for each taxon, using the stdout files
  (3) ./make_plots.sh to make plots (one per taxon) from the summary data
