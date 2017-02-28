'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import datetime
import sys

from Bio.Seq import Seq
from _regex_core import parse_set_member
from synbiochem.utils import seq_utils


_DEFAULT_CODONS = {am_ac: 'NNK' for am_ac in seq_utils.AA_CODES.values()}


def get_oligos(templ, set_len=1, melt_temp=60, codons=None):
    '''Gets oligos.'''
    oligos = []
    def_codons = _get_codons(codons)

    for set_idx, start_pos in enumerate(xrange(0, len(templ), 3 * set_len)):
        oligos.extend(_get_set(templ, set_idx, start_pos, set_len, melt_temp,
                               def_codons))

    return oligos


def _get_codons(codons):
    '''Gets codons.'''
    if codons is None:
        codons = {}

    def_codons = dict(_DEFAULT_CODONS)
    def_codons.update(codons)
    return def_codons


def _get_set(templ, set_idx, start_pos, set_len, melt_temp, def_codons):
    '''Gets a set.'''
    oligos = []

    start_pos = set_idx * set_len * 3
    end_pos = start_pos + 3 * set_len
    pre_seq, pre_tm = _get_seq_by_tm(templ[:start_pos], melt_temp, False)
    post_seq, post_tm = _get_seq_by_tm(templ[end_pos:], melt_temp)

    for set_member in range(set_len):
        oligo = _get_oligo(templ, set_idx, set_member, start_pos,
                           end_pos, pre_seq,
                           post_seq, def_codons)

        oligo.extend([pre_tm, post_tm,
                      str(datetime.date.today()) +
                      '_' + str(set_idx + 1) +
                      '.' + str(set_member + 1)])
        oligos.append(oligo)

    return oligos


def _get_oligo(templ, set_idx, set_member, start_pos, end_pos,
               pre_seq, post_seq, def_codons):
    '''Gets oligo.'''
    # print set_idx, set_member, start_pos, end_pos

    pos = start_pos + 3 * set_member
    codon = Seq(templ[pos:pos + 3])

    seq = pre_seq + \
        templ[start_pos:pos] + \
        def_codons[str(codon.translate())] + \
        templ[pos + 3: end_pos] + \
        post_seq

    return [set_idx + 1,
            set_member + 1,
            seq,
            len(seq)]


def _get_seq_by_tm(seq, melt_temp, forward=True):
    t_m = seq_utils.get_melting_temp(seq) if seq else float('NaN')

    try:
        seq, t_m = seq_utils.get_seq_by_melt_temp(seq, melt_temp, forward)
    except ValueError:
        pass

    return seq, t_m


def main(args):
    '''main method.'''
    for idx, oligo in enumerate(get_oligos(args[0],
                                           set_len=int(args[1]),
                                           melt_temp=float(args[2]),
                                           codons={'S': 'DBK'})):
        print '\t'.join(str(val) for val in [idx + 1] + oligo)


if __name__ == '__main__':
    main(sys.argv[1:])
