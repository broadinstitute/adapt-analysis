"""Make table of virus sequences on GenBank.

The output is meant to plot growth in GenBank viral genomes
over time.

This uses entry creation dates for the year, rather than
sample collection dates, in case there is bias in sample
collection dates (e.g., people did not start using them
until more recent); not every entry has a collection date,
but every entry on GenBank has a entry creation date.

To control for some viruses having many segments (and thus
sequences), this only outputs counts for one segment for
each species; it chooses the segment to be the one with
the most number of sequences.
"""

from collections import Counter
from collections import defaultdict

from adapt.prepare import align
from adapt.prepare import ncbi_neighbors

__author__ = 'Hayden Metsky <hayden@mit.edu>'


TAXONOMIES = 'taxonomies.tsv'
OUTPUT_PER_YEAR = 'data/counts-per-year.tsv'
OUTPUT_CUMULATIVE = 'data/counts-cumulative.tsv'

START_YEAR = 2000
END_YEAR = 2020


def num_unique_kmers(accessions, k=31):
    """Download sequences and count number of unique k-mers.

    This uses k=31 by default, which is what kraken uses for metagenomic
    classification.

    Args:
        accessions: list of accessions
        k: k-mer size

    Returns:
        number of unique k-mers
    """
    accessions = list(accessions)
    seqs_fp = ncbi_neighbors.fetch_fastas(accessions)
    seqs = align.read_unaligned_seqs(seqs_fp)
    all_kmers = set()
    for accver, seq in seqs.items():
        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        all_kmers |= set(kmers)
    seqs_fp.close()

    return len(all_kmers)


def main():
    # Read taxonomies and group by segment
    print("Reading taxonomies")
    taxa = {}
    with open(TAXONOMIES) as f:
        for i, line in enumerate(f):
            if i == 0:
                # Skip header
                continue
            ls = line.split('\t')
            family, genus, species, taxid, segment, refseqs = ls
            taxid = int(taxid)
            if taxid in taxa:
                assert taxa[taxid]['family'] == family
                assert taxa[taxid]['genus'] == genus
                assert taxa[taxid]['species'] == species
            else:
                taxa[taxid] = {'family': family,
                               'genus': genus,
                               'species': species,
                               'segments': []}

            taxa[taxid]['segments'].append(segment)

    # Get list of accessions
    print("Fetching list of accessions")
    taxid_to_skip = []
    n = 0
    for taxid in taxa.keys():
        taxon = taxa[taxid]
        if taxid == 11320 or taxid == 11520:
            # Influenza A or B
            neighbors = ncbi_neighbors.construct_influenza_genome_neighbors(taxid)
        else:
            neighbors = ncbi_neighbors.construct_neighbors(taxid)

        # Determine the segment with the most number of sequences, and
        # use this as the 'representative' segment
        if len(neighbors) == 0:
            # Skip taxid
            taxid_to_skip += [taxid]
        elif taxon['segments'] == ['None']:
            # Not segmented
            taxon['rep_segment'] = None
            taxon['rep_segment_acc'] = [n.acc for n in neighbors]
        else:
            neighbor_segments = [n.segment for n in neighbors]
            c = Counter(neighbor_segments)
            rep_segment = c.most_common(1)[0][0]
            taxon['rep_segment'] = rep_segment
            taxon['rep_segment_acc'] = [n.acc for n in neighbors if n.segment == rep_segment]

        n += 1
        print("  Fetched accessions for {} of {} taxonomies".format(n, len(taxa)))

    for taxid in taxid_to_skip:
        del taxa[taxid]

    # Fetch metadata for the representative segment accessions
    print("Fetching metadata")
    n = 0
    for taxid in taxa.keys():
        taxon = taxa[taxid]
        acc = taxon['rep_segment_acc']
        metadata = ncbi_neighbors.fetch_metadata(acc)
        taxon['metadata'] = metadata

        n += 1
        print("  Fetched metadata for {} for {} taxonomies".format(n, len(taxa)))

    # Fetch sequences and compute number of unique k-mers
    print("Determining number of unique k-mers")
    n = 0
    for taxid in taxa.keys():
        taxon = taxa[taxid]

        # Do this for the collection of accessions from each year
        acc_for_year = defaultdict(set)
        for acc, m in taxon['metadata'].items():
            acc_year = m['entry_create_year']
            acc_for_year[acc_year].add(acc)
        nuk_for_year = defaultdict(int)
        for year, accessions in acc_for_year.items():
            nuk = num_unique_kmers(accessions)
            nuk_for_year[year] = nuk
        taxon['num_unique_kmers'] = nuk_for_year

        n += 1
        print("  Counted unique k-mers for {} for {} taxonomies".format(n, len(taxa)))

    # Write TSV
    print("Writing output TSVs")
    with open(OUTPUT_PER_YEAR, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')

        header = ['taxid', 'family', 'genus', 'species', 'segment', 'year',
                'num_genomes', 'num_unique_kmers']
        write_row(header)

        for taxid in taxa.keys():
            taxon = taxa[taxid]
            year_genome_count = defaultdict(int)
            for acc, m in taxon['metadata'].items():
                year_genome_count[m['entry_create_year']] += 1
            for year in range(START_YEAR, END_YEAR+1):
                row = [taxid, taxon['family'], taxon['genus'],
                        taxon['species'], taxon['rep_segment'],
                        year, year_genome_count[year],
                        taxon['num_unique_kmers'][year]]
                write_row(row)

    with open(OUTPUT_CUMULATIVE, 'w') as fw:
        def write_row(row):
            fw.write('\t'.join(str(x) for x in row) + '\n')

        header = ['taxid', 'family', 'genus', 'species', 'segment', 'year',
                'cumulative_num_genomes', 'cumulative_num_unique_kmers']
        write_row(header)

        for taxid in taxa.keys():
            taxon = taxa[taxid]
            for year in range(START_YEAR, END_YEAR+1):
                year_genome_count = 0
                for acc, m in taxon['metadata'].items():
                    if m['entry_create_year'] <= year:
                        year_genome_count += 1
                year_kmer_count = 0
                for y, nuk in taxon['num_unique_kmers'].items():
                    if y <= year:
                        year_kmer_count += nuk
                row = [taxid, taxon['family'], taxon['genus'],
                        taxon['species'], taxon['rep_segment'],
                        year, year_genome_count, year_kmer_count]
                write_row(row)


if __name__ == "__main__":
    main()
