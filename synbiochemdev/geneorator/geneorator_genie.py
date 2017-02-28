'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys

from Bio.Seq import Seq
from synbiochem.utils import seq_utils

_DEFAULT_CODONS = {am_ac: 'NNK' for am_ac in seq_utils.AA_CODES.values()}


def get_oligos(templ, set_len=1, melt_temp=60, codons=None):
    '''Gets oligos.'''
    oligos = []

    if codons is None:
        codons = {}

    def_codons = dict(_DEFAULT_CODONS)
    def_codons.update(codons)

    for set_idx in xrange(0, len(templ), 3 * set_len):
        pos = set_idx
        codon = Seq(templ[pos:pos + 3])

        pre_seq, pre_tm = _get_seq_by_tm(templ[:pos], melt_temp, False)
        post_seq, post_tm = _get_seq_by_tm(templ[pos + 3:], melt_temp)

        oligos.append([(pos / 3) + 1,
                       pre_seq + def_codons[str(codon.translate())] + post_seq,
                       pre_tm,
                       post_tm])

    return oligos


def _get_seq_by_tm(seq, melt_temp, forward=True):
    t_m = float('NaN')

    try:
        seq, t_m = seq_utils.get_seq_by_melt_temp(seq, melt_temp, forward)
    except ValueError:
        pass

    return seq, t_m


def main(args):
    '''main method.'''
    for oligo in get_oligos(args[0],
                            melt_temp=float(args[1]),
                            codons={'S': 'DBK'}):
        print '\t'.join([str(val) for val in oligo])


if __name__ == '__main__':
    main(sys.argv[1:])
