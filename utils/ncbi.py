"""Utilities for downloading from NCBI.

This is just a wrapper around the ncbi_neighbors module from adapt.
"""

from adapt.prepare import ncbi_neighbors

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def fetch_neighbors_acc(taxid, segment, influenza=False):
    """Fetch sequence accessions from NCBI.

    Args:
        taxid: taxonomic ID to download neighbors for
        segment: only keep sequences labeled with this segment (or
            None if unsegmented)
        influenza: if True, assume taxid represents an influenza
            A or B virus taxonomy, and fetch sequences using NCBI's
            influenza database

    Returns:
        collection of accessions
    """
    # Download neighbors for taxid
    if influenza:
        neighbors = ncbi_neighbors.construct_influenza_genome_neighbors(taxid)
    else:
        neighbors = ncbi_neighbors.construct_neighbors(taxid)

    # Filter neighbors by segment
    if segment != None and segment != '':
        neighbors = [n for n in neighbors if n.segment == segment]

    # An accession can appear multiple times in the table (e.g., for
    # multiple RefSeqs); only return each accession once
    acc = set([n.acc for n in neighbors])
    return list(acc)


# Expose ncbi_neighbors's fetch_xml() function
fetch_xml = ncbi_neighbors.fetch_xml

# Expose ncbi_neighbors's parse_genbank_xml_for_source_features() function
parse_genbank_xml_for_source_features = ncbi_neighbors.parse_genbank_xml_for_source_features

# Expose ncbi_neighbors's fetch_fastas() function
fetch_fastas = ncbi_neighbors.fetch_fastas

# Expose ncbi_neighbors's fetch_metadata() function
fetch_metadata = ncbi_neighbors.fetch_metadata
