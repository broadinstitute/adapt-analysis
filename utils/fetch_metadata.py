"""Download metadata for a list of accessions.
"""

import argparse

import ncbi


def add_species(metadata):
    accessions = metadata.keys()
    species = ncbi.get_subtaxa_groups(accessions, 'species')

    # species is {species: [accessions]}; invert to be {accession: species}
    species_inv = {}
    for s, al in species.items():
        for a in al:
            a = a.split('.')[0] # remove version
            species_inv[a] = s

    for a in accessions:
        metadata[a]['species'] = species_inv[a]
    return metadata


def make_table(metadata, out_tsv):
    with open(out_tsv, 'w') as fw:
        header = ['accession', 'country', 'year', 'entry_create_year',
                'taxid', 'species']
        fw.write('\t'.join(header) + '\n')
        for accession in metadata:
            row = [accession] + [str(metadata[accession][c]) for c in header[1:]]
            fw.write('\t'.join(row) + '\n')


def main(args):
    metadata = ncbi.fetch_metadata(args.a)
    metadata = add_species(metadata)
    make_table(metadata, args.o)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', nargs='+',
            help=("Space-separated list of accessions"))
    parser.add_argument('-o', help=("Output TSV"))
    args = parser.parse_args()

    main(args)
