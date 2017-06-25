'''
Created on 30 Mar 2017

@author: neilswainston
'''
import sys

from synbiochem.utils import seq_utils


def _get_subjects(filename):
    '''Gets subjects.'''
    id_seqs = {}

    with open(filename) as fle:
        for line in fle.read().splitlines():
            tokens = line.split('\t')
            id_seqs[tokens[1]] = tokens[2].upper()

    return id_seqs


def _get_queries(filename):
    '''Gets queries.'''
    id_seqs = {}

    with open(filename) as fle:
        for line in fle.read().splitlines():
            tokens = line.split('\t')
            id_seqs[tokens[1] + ' ' + tokens[3]] = tokens[4].upper()

    return id_seqs


def main(args):
    '''main method.'''
    id_seqs_subjects = _get_subjects(args[0])
    id_seqs_queries = _get_queries(args[1])

    for result in seq_utils.do_blast(id_seqs_subjects, id_seqs_queries):
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                print result.query + '\t' + alignment.title + '\t' + \
                    str(hsp.align_length) + '\t' + str(hsp.gaps) + '\t' + \
                    str(result.query_length) + '\t' + str(alignment.length)
                # print hsp
                # print


if __name__ == '__main__':
    main(sys.argv[1:])
