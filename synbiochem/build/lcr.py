'''
lcr (c) University of Manchester 2015

lcr is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import collections
import itertools
import sys
from synbiochem.build import melting_temp_utils

reagent_concs = melting_temp_utils.DEFAULT_REAGENT_CONCS
melting_temp_calculator = \
    melting_temp_utils.MeltingTempCalculator(reagent_concs)


def get_bridging_oligos(target_melt_temp, sequences, plasmid_seq=None,
                        shuffle=False):
    '''Designs dominos (bridging oligos) for LCR.'''
    num_sequences = len(sequences)
    perms = itertools.permutations(range(num_sequences))
    orderings = sorted([list(a) for a in set(perms)]) \
        if shuffle \
        else [range(num_sequences)]

    if plasmid_seq is not None:
        sequences.append(plasmid_seq)
        for ordering in orderings:
            ordering.insert(0, num_sequences)
            ordering.append(num_sequences)

    print orderings

    pairs = []

    for ordering in orderings:
        pairs.extend(_pairwise(ordering))

    print pairs
    print collections.Counter(pairs)
    print [a for a in sorted(set(pairs))]

    for pair in sorted(set(pairs)):
        reverse, reverse_tm = \
            _get_bridge(sequences[pair[0]], False, target_melt_temp)
        forward, forward_tm = \
            _get_bridge(sequences[pair[1]], True, target_melt_temp)
        print pair[0], pair[1], reverse, forward, reverse_tm, forward_tm, \
            reverse + forward if reverse is not None and forward is not None \
            else ''


def _pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    first, second = itertools.tee(iterable)
    next(second, None)
    return itertools.izip(first, second)


def _get_bridge(sequence, forward, target_melt_temp):
    '''Gets half of bridging oligo.'''
    for i in range(len(sequence)):
        subsequence = sequence[:(i+1)] if forward else sequence[-(i+1):]
        melting_temp = melting_temp_calculator.get_melting_temp(subsequence)

        if melting_temp > target_melt_temp:
            return subsequence, melting_temp

    return None, -float('inf')


def _read_file(filename):
    melting_temp = 70
    sequences = []
    plasmid_seq = None
    shuffle = False

    in_plasmid = False
    seq = None

    with open(filename) as fle:
        for line in fle:
            if line.startswith('MELTING_TEMP:'):
                melting_temp = float(line.replace('MELTING_TEMP:', '').strip())
            elif line.startswith('SHUFFLE:'):
                shuffle = line.replace('MELTING_TEMP:', '').strip() == 'True'
            elif line.startswith('>PLASMID'):
                in_plasmid = True

                if seq is not None:
                    sequences.append(seq)

                seq = []
            elif line.startswith('>'):
                if seq is not None:
                    sequences.append(seq)

                seq = []
            elif len(line) > 0 and seq is not None:
                seq.append(line)

    return melting_temp, sequences, plasmid_seq, shuffle


def main(argv):
    '''main method'''
    melting_temp, sequences, plasmid_seq, shuffle = _read_file(argv[1])
    get_bridging_oligos(melting_temp, sequences, plasmid_seq, shuffle)


if __name__ == '__main__':
    main(sys.argv)

# get_bridging_oligos(70, ['saaaaae', 'sccccce', 'sggggge'], 'sttttte', True)
