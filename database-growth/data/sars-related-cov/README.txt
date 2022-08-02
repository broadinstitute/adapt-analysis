Handled manually by:
- Downloading table for all sequences released up to 12/31/2019 at https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome-related%20coronavirus,%20taxid:694009&CreateDate_dt=1900-01-01T00:00:00.00Z%20TO%202019-12-31T23:59:59.00Z
- Counting number of sequences per year for those sequences (table downloaded in csv format, containing release years)
- Using the web interface to count number of sequences for 2020, 2021, and 2022

Number of SARS-related CoV sequences through 2019, counted from downloaded CSV file:
> cat Downloads/sequences.csv | tail -n +2 | awk -F',' '{print $2}' | sort | awk -F'-' '{print $1}' | sort -n | uniq -c
  80 2003
 344 2004
 106 2005
 132 2006
 127 2007
 132 2008
  64 2009
  81 2010
  52 2011
  36 2012
  76 2013
  42 2014
   8 2015
   2 2016
  27 2017
 102 2018

For 2020: 47290
For 2021: 2977018
For 2022: 2986793
