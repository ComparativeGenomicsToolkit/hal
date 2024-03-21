#!/usr/bin/env python3

"""
Make a table out of halCoverage output

"""

import os, sys
from argparse import ArgumentParser

def main():
    parser = ArgumentParser()

    parser.add_argument('--input', nargs = '+', type=str, required=True, help = 'Input files, must be output of halCoverage')
    parser.add_argument('--reference', type=str, required=True, help = 'Name of reference genome used in halCoverage')
    parser.add_argument('--chroms', nargs = '+', type=str, help = 'Lump all chroms together except these when doing breakdown')
    parser.add_argument('--counts', action='store_true', help = 'Write counts instead of percentages')

    options = parser.parse_args()

    # map filename to coverage tables
    file_coverage = {}

    for input_path in options.input:
        cov_type = None
        name = os.path.splitext(os.path.basename(input_path))[0]
        assert name not in file_coverage
        # map coverage type to table (where type is total or chromosome)
        file_coverage[name] = {}
        with open(input_path, 'r') as input_file:
            for line in input_file:                
                if line.startswith('Genome, sites'):
                    cov_type = 'Total'
                    assert cov_type not in file_coverage[name]
                    file_coverage[name][cov_type] = {}
                elif line.startswith('Coverage on'):
                    cov_type = line.rstrip()[len('Coverage on '):]
                    assert cov_type not in file_coverage[name]
                    file_coverage[name][cov_type] = {}
                    
                else:
                    assert cov_type
                    toks = line.rstrip().replace(' ', '').split(',')
                    if len(toks) > 1:
                        species = toks[0]
                        coverage = int(toks[1])
                        file_coverage[name][cov_type][species] = coverage


        # merge up the chromosomes
        if options.chroms:
            assert 'Chroms' not in file_coverage[name]
            file_coverage[name]['Chroms'] = {}
            to_remove = []
            for chrom_name, chrom_cov in file_coverage[name].items():
                if chrom_name not in ['Chroms', 'Total'] + options.chroms:
                    for species_name, species_cov in chrom_cov.items():
                        if species_name in file_coverage[name]['Chroms']:
                            file_coverage[name]['Chroms'][species_name] += species_cov
                        else:
                            file_coverage[name]['Chroms'][species_name] = species_cov
                    to_remove.append(chrom_name)
            for chrom_name in to_remove:
                del file_coverage[name][chrom_name]
            
                
    # get the labels
    table_names = set(file_coverage.keys())
    chrom_names = set()
    species_names = set()
    for file_name, file_cov in file_coverage.items():
        chroms = set(file_cov.keys())
        if not chrom_names:
            chrom_names = chroms
        else:
            chrom_name = chrom_names.intersection(chroms)
        for chrom_name, chrom_cov in file_cov.items():
            species = set(chrom_cov.keys())
            if not species_names:
                species_names = species
            else:
                species_names = species_names.intersection(species)

    if not options.counts:
        assert options.reference in species_names

    # push totals to front of chroms
    chrom_names_sorted = []
    for chrom in ['Total', 'Chroms']:
        if chrom in chrom_names:
            chrom_names_sorted.append(chrom)
    for chrom in sorted(list(chrom_names)):
        if chrom not in ['Total', 'Chroms']:
            chrom_names_sorted.append(chrom)

    # don't bother with reference species since it'll be 100
    if not options.counts and options.reference:
        species_names.remove(options.reference)
        
    # print the tables
    for chrom_name in chrom_names_sorted:
        sys.stdout.write('Coverage for {}\n'.format(chrom_name))
        # header
        sys.stdout.write('Name')        
        for species in sorted(list(species_names)):
            sys.stdout.write('\t{}'.format(species))
        sys.stdout.write('\n')
        # rows
        for name in sorted(list(table_names)):
            sys.stdout.write(name)
            if not options.counts:
                ref = file_coverage[name][chrom_name][options.reference]
            for species in sorted(list(species_names)):
                cov = file_coverage[name][chrom_name][species]
                if not options.counts:
                    cov = round(100 * (float(cov) / ref), 2)
                sys.stdout.write('\t{}'.format(cov))
            sys.stdout.write('\n')
    
    return 0

if __name__ == "__main__":
    main()
