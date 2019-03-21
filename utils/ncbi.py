"""Utilities for downloading from NCBI.
"""

from collections import defaultdict
import tempfile
import time
import urllib.parse
import urllib.request
from xml.dom import minidom

__author__ = 'Hayden Metsky <hayden@mit.edu>'


def ncbi_neighbors_url(taxid):
    """Construct URL for downloading list of genome neighbors.

    Args:
        taxid: taxonomic ID to download neighbors for

    Returns:
        str representing download URL
    """
    params = urllib.parse.urlencode({'taxid': taxid, 'cmd': 'download2'})
    url = 'https://www.ncbi.nlm.nih.gov/genomes/GenomesGroup.cgi?%s' % params
    return url


def fetch_neighbors_acc(taxid, segment_to_keep=None):
    """Fetch sequence accessions from NCBI.

    Args:
        taxid: taxonomic ID to download neighbors for
        segment_to_keep: only keep sequences labeled with this segment (or
            None if unsegmented)

    Returns:
        collection of accessions
    """
    accs = []
    url = ncbi_neighbors_url(taxid)
    r = urllib.request.urlopen(url)
    raw_data = r.read()
    for line in raw_data.decode('utf-8').split('\n'):
        line_rstrip = line.rstrip()
        if line_rstrip.startswith('## Columns'):
            # Skip header
            continue
        if line_rstrip != '':
            line_rstrip_s = line_rstrip.split('\t')
            if len(line_rstrip_s) > 1:
                acc = line_rstrip_s[1]
                segment = line_rstrip_s[5].replace('segment', '').strip()
                if segment_to_keep is not None and segment_to_keep != segment:
                    # Wrong segment
                    continue
                accs += [acc]
    return accs


def ncbi_fasta_download_url(accessions):
    """Construct URL for downloading FASTA sequence.

    Args:
        accessions: collection of accessions to download sequences for

    Returns:
        str representing download URL
    """
    ids = ','.join(accessions)
    # Use safe=',' to not encode ',' as '%2'
    params = urllib.parse.urlencode({'id': ids, 'db': 'nuccore',
        'rettype': 'fasta', 'retmode': 'text'}, safe=',')
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?%s' % params
    return url


def fetch_fastas(accessions, batch_size=100, reqs_per_sec=2):
    """Download sequences for accessions.

    Entrez enforces a limit of ~3 requests per second (or else it
    will return a 'Too many requests' error); to avoid this, this
    aims for ~2 requests per second. To use up to 10 requests per second,
    request an API key from Entrez.

    Args:
        accessions: collection of accessions to download sequences for
        batch_size: number of accessions to download in each batch
        reqs_per_sec: number of requests per second to allow

    Returns:
        tempfile object containing the sequences in fasta format
    """
    # Make temp file
    fp = tempfile.NamedTemporaryFile()

    # Download sequences in batches
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:(i + batch_size)]
        url = ncbi_fasta_download_url(batch)
        r = urllib.request.urlopen(url)
        raw_data = r.read()
        for line in raw_data.decode('utf-8').split('\n'):
            fp.write((line + '\n').encode())
        time.sleep(1.0/reqs_per_sec)

    # Set position to 0 so it can be re-read
    fp.seek(0)

    return fp


def ncbi_xml_download_url(accessions):
    """Construct URL for downloading GenBank XML.

    Args:
        accessions: collection of accessions to download XML for

    Returns:
        str representing download URL
    """
    ids = ','.join(accessions)
    # Use safe=',' to not encode ',' as '%2'
    params = urllib.parse.urlencode({'id': ids, 'db': 'nuccore',
        'retmode': 'xml'}, safe=',')
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?%s' % params
    return url


def fetch_xml(accessions, batch_size=100, reqs_per_sec=2):
    """Download XML for accessions.

    Entrez enforces a limit of ~3 requests per second (or else it
    will return a 'Too many requests' error); to avoid this, this
    aims for ~2 requests per second. To use up to 10 requests per second,
    request an API key from Entrez.

    Args:
        accessions: collection of accessions to download XML for
        batch_size: number of accessions to download in each batch
        reqs_per_sec: number of requests per second to allow

    Returns:
        tempfile object containing the downloaded XML data
    """
    # Make temp file
    fp = tempfile.NamedTemporaryFile()

    # Only write the header once; otherwise, it will be written for each
    # beach, and then the file will not be able to be parsed
    def is_xml_header(line):
        return (line.startswith('<?xml ') or line.startswith('<!DOCTYPE '))

    # Download in batches
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:(i + batch_size)]
        url = ncbi_xml_download_url(batch)
        r = urllib.request.urlopen(url)
        raw_data = r.read()
        for line in raw_data.decode('utf-8').split('\n'):
            if i > 0 and is_xml_header(line):
                # Only write XML header for the first batch (i == 0)
                continue
            if i > 0 and '<GBSet>' in line:
                # Do not write GBSet open after the first batch
                line = line.replace('<GBSet>', '')
            if '</GBSet>' in line:
                # Never write GBSet close until the end
                line = line.replace('</GBSet>', '')
            fp.write((line + '\n').encode())
        time.sleep(1.0/reqs_per_sec)

    # Write the GBSet close
    fp.write(('</GBSet>' + '\n').encode())

    # Set position to 0 so it can be re-read
    fp.seek(0)

    return fp


def parse_genbank_xml_for_source_features(fn):
    """Parse GenBank XML to extract source features.

    Args:
        fn: path to XML file, as generated by GenBank

    Returns:
        dict {accession: [(qualifier name, qualifier value)]}
    """
    doc = minidom.parse(fn)
    def parse_xml_node_value(element, tag_name):
        return element.getElementsByTagName(tag_name)[0].firstChild.nodeValue

    source_features = defaultdict(list)

    seqs = doc.getElementsByTagName('GBSeq')
    for seq in seqs:
        accession = parse_xml_node_value(seq, 'GBSeq_primary-accession')
        feature_table = seq.getElementsByTagName('GBSeq_feature-table')[0]
        for feature in feature_table.getElementsByTagName('GBFeature'):
            feature_key = parse_xml_node_value(feature, 'GBFeature_key')
            if feature_key == 'source':
                quals = feature.getElementsByTagName('GBFeature_quals')[0]
                for qualifier in quals.getElementsByTagName('GBQualifier'):
                    qual_name = parse_xml_node_value(qualifier, 'GBQualifier_name')
                    qual_value = parse_xml_node_value(qualifier, 'GBQualifier_value')
                    source_features[accession].append((qual_name, qual_value))

    return source_features
