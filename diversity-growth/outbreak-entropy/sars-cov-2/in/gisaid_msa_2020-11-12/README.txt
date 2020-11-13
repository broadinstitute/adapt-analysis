This full alignment (msa_1112.fasta) is based on 196,293 submissions to EpiCoV.
Both duplicate and low-quality sequences (>5% NNNNs) have been removed, using only complete sequences (length >29,000 bp). The resulting alignment of 184,197 sequences is created using mafft https://doi.org/10.1093/molbev/mst010 in 3 separate steps.

(1) Each sequence is individually aligned to the reference hCoV-19/Wuhan/WIV04/2019 (EPI_ISL_402124). Sequences that created dubious insertions of >=30 nucleotides in the reference sequence and occurred only once in the database are discarded. The alignments are created with the command:
mafft --thread -1 input.fasta > output.fasta

(2) All sequences that result in insertions in the reference from (1) are then aligned with a opening gap penalty of 10 to prevent long stretches of dubious insertions in the alignment due to the presence of long stretches of NNNNs. The following command is used:
mafft --retree 3 --maxiterate 10 --thread -1 --nomemsave --op 10 seqsCausingInsertionsInRef.fasta > seqsCausingInsertionsInRef_aligned.fasta

(3) The rest of the sequences that did not result in insertions are aligned to the resulting alignment in step (2) with this command:
mafft --thread 1 --quiet --keeplength --add sequencesNotCausingInsertionsInRef.fa seqsCausingInsertionsInRef_aligned.fasta > msa_1112.fasta

Thanks to Rob Lanfear for discussion and suggestions.