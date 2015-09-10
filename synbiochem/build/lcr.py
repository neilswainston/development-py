'''
lcr (c) University of Manchester 2015

lcr is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
import sys
from synbiochem.build import melting_temp_utils


def get_bridging_oligos(target_melt_temp, sequences, plasmid_seq=None,
                        shuffle=False):
    '''Designs dominos (bridging oligos) for LCR.'''
    num_sequences = len(sequences)
    orderings = sorted([list(a) for a in \
                        set(itertools.permutations(range(num_sequences)))]) \
        if shuffle \
        else [range(num_sequences)]

    if plasmid_seq is not None:
        sequences.append(plasmid_seq)
        for ordering in orderings:
            ordering.insert(0, num_sequences)
            ordering.append(num_sequences)

    pairs = []

    for ordering in orderings:
        pairs.extend(_pairwise(ordering))

    melting_temp_calc = melting_temp_utils.MeltingTempCalculator()

    result_file = open('result.xls', 'w+')

    for pair in sorted(set(pairs)):
        reverse, reverse_tm = _get_bridge(sequences[pair[0]], False,
                                          target_melt_temp, melting_temp_calc)
        forward, forward_tm = _get_bridge(sequences[pair[1]], True,
                                          target_melt_temp, melting_temp_calc)
        result = [str(pair[0]) if pair[0] is not num_sequences else 'P', \
            str(pair[1]) if pair[1] is not num_sequences else 'P', \
            reverse, forward, str(reverse_tm), str(forward_tm), \
            reverse + forward if reverse is not None and forward is not None \
            else '']

        result_line = '\t'.join(result)
        print result_line
        result_file.write(result_line + '\n')

    result_file.close()


def _pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    first, second = itertools.tee(iterable)
    next(second, None)
    return itertools.izip(first, second)


def _get_bridge(sequence, forward, target_melt_temp, melting_temp_calculator):
    '''Gets half of bridging oligo.'''
    for i in range(len(sequence)):
        subsequence = sequence[:(i + 1)] if forward else sequence[-(i + 1):]
        melting_temp = melting_temp_calculator.get_melting_temp(subsequence)

        if melting_temp > target_melt_temp:
            return subsequence, melting_temp

    return None, -float('inf')


def _read_parameter_file(filename):
    '''Reads parameter file.'''
    melting_temp = 70
    sequences = []
    plasmid_seq = None
    shuffle = False

    in_plasmid = False
    seq = ''

    with open(filename) as fle:
        for line in fle.read().splitlines():
            if line.startswith('MELTING_TEMP:'):
                melting_temp = float(line.replace('MELTING_TEMP:', '').strip())
            elif line.startswith('SHUFFLE:'):
                shuffle = line.replace('MELTING_TEMP:', '').strip() == 'True'
            elif line.startswith('>'):
                if len(seq) > 0:
                    if in_plasmid:
                        plasmid_seq = seq
                    else:
                        sequences.append(seq)

                    seq = ''

                if line.startswith('>PLASMID'):
                    in_plasmid = True
            elif len(line) > 0 and not line.startswith('#') and seq is not None:
                seq += line.upper()

        if seq is not None:
            if in_plasmid:
                plasmid_seq = seq
                in_plasmid = False
            else:
                sequences.append(seq)

    return melting_temp, sequences, plasmid_seq, shuffle


def main(argv):
    '''main method'''
    melting_temp, sequences, plasmid_seq, shuffle = \
        _read_parameter_file(argv[1])
    get_bridging_oligos(melting_temp, sequences, plasmid_seq, shuffle)


if __name__ == '__main__':
    main(sys.argv)

# get_bridging_oligos(70, ['saaaaae', 'sccccce', 'sggggge'], 'sttttte', True)
