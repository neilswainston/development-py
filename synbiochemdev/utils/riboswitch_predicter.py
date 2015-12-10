'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
from synbiochem.utils.sequence_utils import CodonOptimiser, \
    get_minimum_free_energy


def main(argv):
    '''main method'''
    upstream_seq = argv[1]
    upstream_trunc_seq = upstream_seq[-int(argv[2]):]
    variant_seq = argv[3]
    downstream_seq = argv[4]

    cod_opt = CodonOptimiser('9606')
    sequences = []

    for rev_trans in cod_opt.get_all_rev_trans(variant_seq):
        sequences.extend([upstream_seq + rev_trans + downstream_seq,
                          upstream_trunc_seq + rev_trans + downstream_seq])

    mfes = get_minimum_free_energy(sequences)

    outfile = open(argv[5], 'w')

    for i in xrange(0, len(sequences), 2):
        outfile.write('\t'.join([sequences[i], sequences[i + 1],
                                 str(mfes[i]), str(mfes[i + 1]),
                                 str(mfes[i] - mfes[i + 1])]) + '\n')

    outfile.close()

if __name__ == '__main__':
    main(sys.argv)
