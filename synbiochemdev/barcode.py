'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import sys

from synbiochem.utils import seq_utils


def get_seqs(num_barcodes, length, max_repeat_nuc):
    '''Get barcodes.'''
    seqs = [seq_utils.get_random_dna(length, max_repeat_nuc)
            for _ in range(num_barcodes)]

    for query_idx, seq in enumerate(seqs[:-1]):
        id_seqs_queries = {query_idx: seq}
        id_seqs_subjects = {'seq_' + str(idx + query_idx + 1): seq
                            for idx, seq in enumerate(seqs[query_idx + 1:])}

        do_blast(id_seqs_subjects, id_seqs_queries)

    return seqs


def do_blast(subjects, queries, exp_threshold=1):
    '''Runs BLAST, filtering results.'''
    for result in seq_utils.do_blast(subjects, queries, word_size=4):
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < exp_threshold:
                    print result.query + '\t' + alignment.hit_def + '\n' + str(hsp)
                    print
                    return False

    return True


def main(args):
    '''main method.'''
    get_seqs(int(args[0]), int(args[1]), int(args[2]))


if __name__ == '__main__':
    main(sys.argv[1:])
