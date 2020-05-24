There are two related analyses here:
  * Entropy: one is meant to show the increase in entropy (nucleotide) along the genome
    of k-mers within a sliding window. If we were to look within a window and choose
    a k-mer (say, 28-mer) to use as a guide for detection, it's to show that the
    diversity of them is growing throughout the genome over time -- and thus we have
    to be careful in how we select it.
  * Variants: the other is meant to show variants along the genome more directly,
    at different points in time. It takes a collection of genomes and 'collapses'
    them down to one, showing a dot at a site if there are genomes with a variant
    at that site. It is to show that the number of variants throughout the genome
    increases over time.

SARS-CoV-2 is a good example because there are so many genomes over the course
of the outbreak, so we can show this as a detailed time-resolution.
