"""Create table of species to download and download FASTAs of genomes
for these species.
"""

import argparse
import datetime
from collections import defaultdict
import os
import re
import statistics
import time

from Bio import Entrez

__author__ = 'Hayden Metsky <hayden@mit.edu>'


Entrez.email = "hayden@mit.edu"


class SequenceFromAccessionList:

    def __init__(self, representative, neighbor, host,
                 lineage, taxonomy_name, segment):
        self.representative = representative
        self.neighbor = neighbor
        self.taxonomy_name = taxonomy_name

        # Parse hosts
        self.hosts = host.lower().split(',')
        self.human_is_host = 'human' in self.hosts

        # Parse lineage
        lineage_split = lineage.split(',')
        num_ranks = len(lineage_split)
        assert num_ranks <= 3
        if num_ranks == 1:
            # Assume only species is given
            self.lineage = ('', '', lineage_split[0])
        elif num_ranks == 2:
            # Assume only genus and species is given
            self.lineage_genus, self.lineage_species = lineage_split
            self.lineage = ('', lineage_split[0], lineage_split[1])
        else:
            # family, genus, and species are given
            self.lineage = (lineage_split[0], lineage_split[1], lineage_split[2])

        # Parse segment
        self.is_segmented = segment != 'segment'
        if not self.is_segmented:
            self.segment = ''
        else:
            assert segment.startswith('segment ')
            assert segment != 'segment '
            self.segment = segment.replace('segment ', '')

        # Group by (lineage, segment)
        self.group = (self.lineage, self.segment)

    def __hash__(self):
        return hash(self.neighbor)

    def __eq__(self, other):
        return self.neighbor == other.neighbor

    @staticmethod
    def from_line(line):
        ls = line.split('\t')

        segment = ls[5]
        removals = ['RNA', 'DNA']
        for removal in removals:
            if (segment.startswith('segment ' + removal + ' ') and
                    len(segment) > len('segment ' + removal + ' ')):
                # change 'segment removal X' to 'segment X'
                segment = 'segment ' + segment[len('segment ' + removal + ' '):]

        return SequenceFromAccessionList(ls[0], ls[1], ls[2], ls[3], ls[4],
                                         segment)


def read_genome_accession_list(fn):
    sequences = []
    with open(fn) as f:
        for line in f:
            if line.startswith('#'):
                # Skip comments
                continue
            sequences += [SequenceFromAccessionList.from_line(line.rstrip())]
    return sequences


def uniqueify_genome_accession_list(sequences):
    # Some sequences have the same name (neighbor) but different representatives
    # (i.e., different reference sequences) and therefore appear more than once
    # in the list; the sequences of these genomes are all we care about, so
    # "unique-ify" them so that each neighbor appears just once
    # (arbitrarily select the sequence, so in effect arbitrarily pick the
    # representative)
    return list(set(sequences))


def filter_sequences_with_nonhuman_host(sequences, args):
    # Return only those sequences that are from a lineage that has
    # at least one sequence for which human is a host

    # Find all lineages that have a sequence for which human is a host
    human_host_lineages = set(s.lineage for s in sequences if s.human_is_host)

    # If provided, add in lineages listed in a file
    if args.human_host_lineages_to_add:
        with open(args.human_host_lineages_to_add) as f:
            for line in f:
                ls = line.rstrip().split('\t')
                human_host_lineages.add(tuple(ls))

    return [s for s in sequences if s.lineage in human_host_lineages]


def download_raw_from_genbank(sequences,
                              results_type='fasta',
                              batch_size=50,
                              max_tries=5,
                              base_delay=5):
    # Entrez gives sporadic exceptions (usually RuntimeErrors);
    # retry the call a few times if this happens before crashing
    # the entire program
    try_num = 1
    while try_num <= max_tries:
        try:
            return _download_raw_from_genbank(sequences,
                                              results_type=results_type,
                                              batch_size=batch_size)
        except Exception as e:
            if try_num == max_tries:
                # used up all tries
                raise e
            time.sleep(2**(try_num - 1) * base_delay)
            try_num += 1


def _download_raw_from_genbank(sequences,
                               results_type='fasta',
                               batch_size=50):
    # Download raw data from GenBank of type 'gb' or 'fasta', as specified
    # by results_type
    accessions = [s.neighbor for s in sequences]

    # Query GenBank using the accessions, fetching results in batches
    # of size batch_size
    raw_results = ''
    for batch_start in range(0, len(accessions), batch_size):
        batch_end = min(batch_start + batch_size, len(accessions))
        batch_accessions = accessions[batch_start:batch_end]

        query = ','.join(batch_accessions)
        reader = Entrez.read(Entrez.epost(db='nuccore', id=query))
        results = Entrez.efetch(db='nuccore',
                                rettype=results_type,
                                retmax=batch_size,
                                webenv=reader['WebEnv'],
                                query_key=reader['QueryKey'])
        raw_results += results.read()

    return raw_results


def parse_metadata_from_gb_results(gb_results):
    # gb_results is output of download_raw_from_genbank with results_type='gb'
    # returns a dict {accession: {metadata_key: metadata_value}}

    # Each result is separated by a line with only '//', so split on this
    gb_results_split = re.split('^//$', gb_results, flags=re.MULTILINE)

    # Setup regex patterns for gb file
    patterns = {'accession': re.compile('^ACCESSION\s+(\w+)( |$)', re.MULTILINE),
                'strain': re.compile('/strain="(.+?)"', re.MULTILINE),
                'isolate': re.compile('/isolate="(.+?)"', re.MULTILINE),
                'collection_date': re.compile('/collection_date="(.+?)"', re.MULTILINE),
                'seq_length': re.compile('^LOCUS\s+\w+\s+(\d+) bp', re.MULTILINE)}

    # Assume year is the 4-digit number in collection_date
    year_pattern = re.compile('([0-9]{4})')

    metadata = {}
    for result in gb_results_split:
        if len(result) == 0 or result.isspace():
            continue

        # Find accession number
        accession_match = patterns['accession'].search(result)
        if not accession_match:
            raise Exception("Unknown accession number in result")
        accession = accession_match.group(1)

        # Find metadata
        result_metadata = {}
        for k in ['strain', 'isolate', 'collection_date', 'seq_length']:
            match = patterns[k].search(result)
            if match:
                result_metadata[k] = match.group(1)
            else:
                result_metadata[k] = None

        # Parse collection_date for a year
        if result_metadata['collection_date'] is not None:
            year_match = year_pattern.search(result_metadata['collection_date'])
            if year_match:
                result_metadata['year'] = int(year_match.group(1))
            else:
                result_metadata['year'] = None
        else:
            result_metadata['year'] = None

        # Parse sequence length
        assert result_metadata['seq_length'] is not None
        result_metadata['seq_length'] = int(result_metadata['seq_length'])

        # Check that the accession hasn't been already encountered
        assert accession not in metadata

        metadata[accession] = result_metadata
    return metadata


def extract_accession_from_header(header):
    acc_match = re.search(
        '\|(?:gb|emb|dbj|ref)\|(.+?)\||gb:(.+?)\|', header)
    if acc_match is None:
        raise Exception("In '%s', could not determine accession" %
                        header)

    if acc_match.group(1):
        return acc_match.group(1)
    else:
        return acc_match.group(2)


present_year = datetime.datetime.now().year
def bin_sequence_years(sequences, metadata, max_years_ago=500,
                       bins=[(0, 5), (5, 10), (10, 15), (15, 20), (20, 500)]):
    # Bin the sequences' years by how many years ago they were
    # to the present year
    bin_num = {b: 0 for b in bins}
    bin_num['unknown'] = 0
    for s in sequences:
        year = metadata[s.neighbor]['year']
        if year is None:
            bin_num['unknown'] += 1
            continue
        years_ago = present_year - year
        if years_ago < 0 or years_ago >= max_years_ago:
            bin_num['unknown'] += 1
            continue

        in_a_bin = False
        for b in bins:
            if years_ago >= b[0] and years_ago < b[1]:
                bin_num[b] += 1
                in_a_bin = True
                break
        assert in_a_bin

    bin_frac = {k: bin_num[k]/float(len(sequences)) for k in bin_num.keys()}
    return bin_frac


def make_species_list(args):
    # Read/parse accession list
    sequences = read_genome_accession_list(args.accession_list)
    sequences = filter_sequences_with_nonhuman_host(sequences, args)
    sequences = uniqueify_genome_accession_list(sequences)

    # Fetch the year for all accessions
    dl = download_raw_from_genbank(sequences, results_type='gb')
    metadata = parse_metadata_from_gb_results(dl)

    # Group sequences by their group field
    by_group = defaultdict(set)
    for s in sequences:
        by_group[s.group].add(s)

    with open(args.output, 'w') as f:
        # Write the header
        header = ['family', 'genus', 'species', 'segment', 'refseqs',
                  'num_sequences', 'seq_length_mean', 'seq_length_cv',
                  'frac_<5yrs', 'frac_5-10yrs', 'frac_10-15yrs',
                  'frac_15-20yrs', 'frac_>20yrs', 'frac_unknown',
                  'sequence_names']
        f.write('\t'.join(header) + '\n')

        for group in sorted(by_group.keys()):
            lineage, segment = group
            if segment is None or segment == '':
                segment = 'N/A'
            family, genus, species = lineage
            group_sequences = by_group[group]
            num_genomes = len(group_sequences)

            # Check lineage fields
            if family == '':
                family = 'unknown'
            if genus == '':
                genus = 'unknown'
            assert species != ''

            # Pull out the RefSeq(s)
            refseqs = set(s.representative for s in group_sequences)
            refseqs_str = ','.join(sorted(refseqs))

            # Pull out the taxonomy names
            taxonomy_names = set(s.taxonomy_name for s in group_sequences)
            taxonomy_names_str = ','.join(sorted(taxonomy_names))

            # Calculate mean and stdev of sequence lengths
            seq_lengths = [metadata[s.neighbor]['seq_length'] for s in group_sequences]
            seq_length_mean = statistics.mean(seq_lengths)
            if len(seq_lengths) > 1:
                seq_length_cv = statistics.stdev(seq_lengths) / float(seq_length_mean)
            else:
                seq_length_cv = 0

            # Bin the years by how many years ago they were to the present
            # year
            year_bin_frac = bin_sequence_years(group_sequences, metadata)

            cols = [family, genus, species, segment, refseqs_str,
                    num_genomes, seq_length_mean, seq_length_cv,
                    year_bin_frac[(0, 5)], year_bin_frac[(5, 10)],
                    year_bin_frac[(10, 15)], year_bin_frac[(15, 20)],
                    year_bin_frac[(20, 500)], year_bin_frac['unknown'],
                    taxonomy_names_str]
            f.write('\t'.join(str(x) for x in cols) + '\n')


def make_fasta_files(args):
    # Read/parse accession list
    sequences = read_genome_accession_list(args.accession_list)
    sequences = filter_sequences_with_nonhuman_host(sequences, args)
    sequences = uniqueify_genome_accession_list(sequences)

    # Group sequences by their group field
    by_group = defaultdict(set)
    for s in sequences:
        by_group[s.group].add(s)

    # Download a set of fasta files for each group of sequences
    for group in sorted(by_group.keys()):
        lineage, segment = group
        if segment is None or segment == '':
            segment = 'NA'
        family, genus, species = lineage

        # Construct a filename for the sequences of this group, and
        # verify it does not exist
        r = ((' ', '_'), ('/', '!'))
        for find, replace in r:
            family = family.replace(find, replace)
            genus = genus.replace(find, replace)
            species = species.replace(find, replace)
            segment = segment.replace(find, replace)
        fn = family + '---' + genus + '---' + species + '---' + segment + '.fasta'
        path = os.path.join(args.output, fn)
        if os.path.exists(path):
            raise Exception("fasta file '%s' already exists" % path)

        # Download the sequences for this group
        group_sequences = by_group[group]
        dl = download_raw_from_genbank(group_sequences, results_type='fasta')

        # Write the fasta
        with open(path, 'w') as f:
            f.write(dl)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    # 'make-species-list' command
    parser_msl = subparsers.add_parser('make-species-list',
        help="Make a list of species")
    parser_msl.add_argument('accession_list',
        help="Viral accession list")
    parser_msl.add_argument('-o', '--output', required=True,
        help="Output TSV of species (split by segment)")
    parser_msl.add_argument('--human-host-lineages-to-add',
        help=("File listing lineages to explicitly include as having "
              "human as a host; each row gives a lineage, tab-separated "
              "by family/genus/species"))
    parser_msl.set_defaults(func=make_species_list)

    # 'make-fasta-files' command
    parser_mff = subparsers.add_parser('make-fasta-files',
        help="Make fasta files of sequences, grouped by lineage and segment")
    parser_mff.add_argument('accession_list',
        help="Viral accession list")
    parser_mff.add_argument('--human-host-lineages-to-add',
        help=("File listing lineages to explicitly include as having "
              "human as a host; each row gives a lineage, tab-separated "
              "by family/genus/species"))
    parser_mff.add_argument('-o', '--output', required=True,
        help="Output directory in which to place fasta files")
    parser_mff.set_defaults(func=make_fasta_files)

    args = parser.parse_args()  
    args.func(args)
