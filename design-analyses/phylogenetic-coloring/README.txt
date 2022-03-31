Construct phylogeny from the MSA used/output by ADAPT, and display this alongside
per-sequence design performance information. This could visualize, for example,
if there is a considerable clade or outlier for this designs are predicted to
not detect the sequence(s) well.

Implementation notes:
* To make tree, I used iqtree with `-m TEST` for model selection
* Tree rooting
    * For LASV, following PMC4537774, I rooted the tree on the 1969 _Pinneo_ strain (KM822127). I did this in FigTree and exported the .tree file in Newick format.
* Plotting
    * Custom conda environment: `conda activate misc-bioinformatic-tools`

