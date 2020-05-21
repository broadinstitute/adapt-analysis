# Uses: https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer

from dna_features_viewer import BiopythonTranslator

# Read GenBank file
graphic_record = BiopythonTranslator().translate_record("NC_045512.2.gb")

# Plot and manually set ticks
ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
ax.set_xticks([0,10000,20000])
ax.figure.savefig("genome.pdf")
