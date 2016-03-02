'''
Created on 1 Mar 2016

@author: neilswainston
'''
import os
import sys
import synbiochem.utils.sequence_utils as sequence_utils
import synbiochemdev.build.lcr as lcr


def main(argv):
    id_seq = {}
    all_oligos = []

    with open(argv[0]) as partsfile:
        next(partsfile)

        for line in partsfile:
            line = line.split('\t')
            id_seq[line[0].strip()] = line[7].upper().strip().replace(
                'gaattcaaaagatctgagtc'.upper(), '').replace(
                'gactcggatccaaactcgag'.upper(), '')

    with open(argv[1]) as designfile:
        for line in designfile:
            line = line.split()[1:]
            pairs = [list(pair) for pair in lcr._pairwise(line)]
            oligos = lcr.get_bridging_oligos(70, [id_seq[val] for val in line],
                                             reagent_concs={sequence_utils.MG:
                                                            0.01})
            all_oligos.extend([pair + oligo[2:]
                               for pair, oligo in zip(pairs, oligos)])

    result_file = open(os.path.splitext(argv[1])[0] + '_dominoes.xls', 'w+')

    for oligo in all_oligos:
        result_line = '\t'.join(oligo)
        print result_line
        result_file.write(result_line + '\n')

    result_file.close()

if __name__ == '__main__':
    main(sys.argv[1:])
