'''
Created on 1 Mar 2016

@author: neilswainston
'''
import os
import sys

import synbiochem.utils
import synbiochemdev.build.lcr as lcr


def make_dominoes(parts_filename, design_filename):
    '''Designs dominoes.'''
    id_seq = {}
    pair_oligos = {}

    with open(parts_filename) as partsfile:
        next(partsfile)

        for line in partsfile:
            line = line.split('\t')
            id_seq[line[0].strip()] = line[7].upper().strip().replace(
                'GAATTCAAAAGATCTGAGTCTTGTA', '').replace(
                'TTGTAGACTCGGATCCAAACTCGAG', '')

    with open(design_filename) as designfile:
        for line in designfile:
            line = line.split()[1:]
            line = line + [line[0]]
            pairs = [pair for pair in synbiochem.utils.pairwise(line)]
            oligos = lcr.get_bridging_oligos(70, [id_seq[val] for val in line])
            pair_oligos.update({pair: oligo[2:]
                                for pair, oligo in zip(pairs, oligos)})

    result_file = open(os.path.splitext(design_filename)[0] + '_dominoes.xls',
                       'w+')

    oligo_pairs = {}

    for pair, oligo in pair_oligos.iteritems():
        result_line = '\t'.join(list(pair) + oligo)
        print result_line
        result_file.write(result_line + '\n')

        if oligo[4] not in oligo_pairs:
            oligo_pairs[oligo[4]] = [list(pair)]
        else:
            oligo_pairs[oligo[4]].append(list(pair))

    result_file.write('\n')

    for oligo, pairs in oligo_pairs.iteritems():
        result_file.write('\t'.join([oligo] + [pair[0] + '_' + pair[1]
                                               for pair in pairs]) + '\n')

    result_file.close()


def main(args):
    '''main method.'''
    for design_filename in args[1:]:
        make_dominoes(args[0], design_filename)

if __name__ == '__main__':
    main(sys.argv[1:])
