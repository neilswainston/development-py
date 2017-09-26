'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
from synbiochem.utils import seq_utils

seqs = [seq_utils.get_random_dna(24, 3) for _ in range(500)]

for query_idx, seq in enumerate(seqs[:-1]):
    id_seqs_subjects = {'seq_' + str(idx + query_idx + 1): seq
                        for idx, seq in enumerate(seqs[query_idx + 1:])}

    for result in seq_utils.do_blast(id_seqs_subjects, {query_idx: seq},
                                     word_size=4):
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                print result.query + '\t' + alignment.hit_def + '\n' + str(hsp)
                print
