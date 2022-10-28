Construct phylogeny from the MSA used/output by ADAPT, and display this alongside
per-sequence design performance information. This could visualize, for example,
if there is a considerable clade or outlier for this designs are predicted to
not detect the sequence(s) well.

Implementation notes:
* To make trees, I used iqtree with `-m TEST` for model selection
* Tree rooting
    * For LASV, following PMC4537774, I rooted the tree on the 1969 _Pinneo_ strain (KM822127). I did this in FigTree and exported the .tree file in Newick format.
    * For Filoviridae, following ICTV (https://talk.ictvonline.org/ictv-reports/ictv_online_report/negative-sense-rna-viruses/w/filoviridae), I rooted the tree on Huangjiao virus (MG599981). I did this in FigTree and export the .tree file in Newick format.
    * For DENV, following PMC6731059 and Nextstrain, I rooted the tree on serotype 4 of dengue virus. I did this in FigTree and exported the .tree file in Newick format. The 4 main clusters, from top to bottom, correspond to serotypes 1, 3, 2, and 4 in that order.
* Plotting
    * Custom conda environment: `conda activate misc-bioinformatic-tools`
* For Filoviridae, all designs were made using the --weight-by-log-size-of-subtaxa argument in ADAPT, except `filoviridae--unweighted`, which did not use that argument.
