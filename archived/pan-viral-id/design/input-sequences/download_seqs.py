"""Create table of species to download and download FASTAs of genomes
for these species.
"""

import argparse
import datetime
from collections import defaultdict
from collections import OrderedDict
import os
import re
import statistics
import textwrap
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


def read_influenza_seqs(fn):
    required_cols = ['accession', 'host', 'segment', 'name']
    possible_segments = [str(i) for i in range(1, 9)]

    sequences = []
    with open(fn) as f:
        header_idx = {}
        for i, line in enumerate(f):
            ls = line.rstrip().split('\t')
            if i == 0:
                # Read header, and check that it contains the required
                # columns
                header_idx = {ls[j]: j for j in range(len(ls))}
                for c in required_cols:
                    assert c in header_idx
                continue

            # Parse row
            neighbor = ls[header_idx['accession']]
            host = ls[header_idx['host']]
            segment = ls[header_idx['segment']].split(' ')[0]
            assert segment in possible_segments
            segment = 'segment ' + segment
            name = ls[header_idx['name']]
            if name.startswith('Influenza A virus'):
                lineage = 'Orthomyxoviridae,Alphainfluenzavirus,Influenza A virus'
            elif name.startswith('Influenza B virus'):
                lineage = 'Orthomyxoviridae,Betainfluenzavirus,Influenza B virus'
            elif name.startswith('Influenza C virus'):
                lineage = 'Orthomyxoviridae,Gammainfluenzavirus,Influenza C virus'
            else:
                raise Exception("Unknown lineage for '%s'" % name)

            sequences += [SequenceFromAccessionList('NA', neighbor, host,
                lineage, name, segment)]

    return sequences


def read_accessions_to_reverse_complement(fn):
    acc_to_rc = set()
    with open(fn) as f:
        for line in f:
            acc_to_rc.add(line.rstrip())
    return acc_to_rc


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

    # Check that all accessions are present in the result (sometimes the
    # result is missing one, which leads to issues later on); if one
    # is missing, raise an Exception so this is re-tried
    if results_type == 'fasta':
        accession_pattern = re.compile('^>(\S+)(?: |$)', re.MULTILINE)
        accessions_found = set(accession_pattern.findall(raw_results))
        # Remove version from accession
        accessions_found = set(a.split('.')[0] for a in accessions_found)
        for a in accessions:
            if a not in accessions_found:
                raise Exception("Accession '%s' is not in results" % a)
    elif results_type == 'gb':
        accession_pattern = re.compile('^ACCESSION\s+(\w+)(?: |$)', re.MULTILINE)
        accessions_found = set(accession_pattern.findall(raw_results))
        for a in accessions:
            if a not in accessions_found:
                raise Exception("Accession '%s' is not in results" % a)

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


def reverse_complement_sequences(acc_to_rc, fasta):
    # acc_to_rc is a collection of accessions whose sequence should
    # be reverse complemented; fasta is a string representing a
    # fasta file
    # This returns fasta, with the sequences in acc_to_rc reverse
    # complemented

    rc_map = {'A': 'T',
              'C': 'G',
              'G': 'C',
              'T': 'A',
              'U': 'A',
              'R': 'Y',
              'Y': 'R',
              'S': 'S',
              'W': 'W',
              'K': 'M',
              'M': 'K',
              'B': 'V',
              'D': 'H',
              'H': 'D',
              'V': 'B',
              'N': 'N',
              '-': '-'}
    def rc(s):
        return ''.join(rc_map.get(b, b) for b in s[::-1])

    # Read a dict of the sequences in fasta
    seqs = OrderedDict()
    curr_seq_name = ""
    for line in fasta.split('\n'):
        line = line.rstrip()
        if len(line) == 0:
            continue
        if curr_seq_name == "":
            # Must encounter a new sequence
            assert line.startswith('>')
        if len(line) == 0:
            # Reset the sequence being read on an empty line
            curr_seq_name = ""
        elif line.startswith('>'):
            curr_seq_name = line[1:]
            seqs[curr_seq_name] = ''
        else:
            # Append the sequence
            seqs[curr_seq_name] += line.upper()

    # Find all sequences whose accession matches one in acc_to_rc,
    # and replace these with their reverse complement; also change
    # the sequence name to reflect this
    accession_pattern = re.compile('^(\S+)(?: |$)')
    seqs_rewritten = OrderedDict()
    for seq_name, seq in seqs.items():
        acc = accession_pattern.match(seq_name).group(1)
        if acc in acc_to_rc:
            seq_name = seq_name + ' [reverse-complement]'
            seq = rc(seq)
        seqs_rewritten[seq_name] = seq
    
    # Re-create the fasta string with the new sequences
    chars_per_line = 70
    fasta_rewritten = ''
    for seq_name, seq in seqs_rewritten.items():
        fasta_rewritten += '>' + seq_name + '\n'
        seq_wrapped = textwrap.wrap(seq, chars_per_line)
        for seq_line in seq_wrapped:
            fasta_rewritten += seq_line + '\n'
        fasta_rewritten += '\n' 

    return fasta_rewritten


def make_species_list(args):
    # Read/parse accession list
    sequences = read_genome_accession_list(args.accession_list)
    if args.influenza_seqs:
        sequences += read_influenza_seqs(args.influenza_seqs)
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
    if args.influenza_seqs:
        sequences += read_influenza_seqs(args.influenza_seqs)
    sequences = filter_sequences_with_nonhuman_host(sequences, args)
    sequences = uniqueify_genome_accession_list(sequences)

    # Read list of sequences to take the reverse complement of
    if args.reverse_complement_accessions:
        acc_to_rc = read_accessions_to_reverse_complement(
            args.reverse_complement_accessions)
    else:
        acc_to_rc = None

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

        # Take the reverse complement of specified sequences
        if acc_to_rc is not None:
            dl = reverse_complement_sequences(acc_to_rc, dl)

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
    parser_msl.add_argument('--influenza-seqs',
        help=("TSV file giving Influenza virus accessions to use; these "
              "may come from the NCBI Influenza Virus Database rather than "
              "the viral genome accession list. First row must be header"))
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
    parser_mff.add_argument('--influenza-seqs',
        help=("TSV file giving Influenza virus accessions to use; these "
              "may come from the NCBI Influenza Virus Database rather than "
              "the viral genome accession list. First row must be header"))
    parser_mff.add_argument('--reverse-complement-accessions',
        help=("File listing accessions (one per row) corresponding to "
              "sequences for which the reverse complement should be taken; "
              "the reverse complement of these sequences will be written "
              "instead"))
    parser_mff.add_argument('-o', '--output', required=True,
        help="Output directory in which to place fasta files")
    parser_mff.set_defaults(func=make_fasta_files)

    args = parser.parse_args()  
    args.func(args)
