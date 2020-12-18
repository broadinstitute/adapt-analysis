sars-related-cov_upto_2019-11-30.fasta.gz is the result of all sequences
in the SARS-related CoV species on NCBI prior to Nov. 30, 2019. In
particular, all sequences with taxid:694009 and nucleotide
completeness="complete" and release date up to 2019-11-30. See:
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome-related%20coronavirus,%20taxid:694009&Completeness_s=complete&CreateDate_dt=1900-01-01T00:00:00.00Z%20TO%202019-11-30T23:59:59.00Z

The sequences include bat SARS-like virus ZC45 and ZXC21 (MG772933 and MG772934).
The latest release date of the sequences is 2018-12-31. Hence, we get the same
sequences if we look up to January 1, 2019. See:
https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome-related%20coronavirus,%20taxid:694009&Completeness_s=complete&CreateDate_dt=1900-01-01T00:00:00.00Z%20TO%202019-01-01T23:59:59.00Z

sars-related-cov_upto_2019-11-30.single-sars-cov-1.fasta.gz is the same
as the above, except downsampled to have just one SARS-CoV-1 genome
(the RefSeq, NC_004718). In particular, I removed all sequences (except
that) whose name starts with 'SARS' -- i.e., the FASTA header contains
the string '|SARS'; these correspond to SARS-CoV-1 on a phylogeny.
