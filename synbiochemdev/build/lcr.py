'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
import os
import sys

import synbiochem.utils
import synbiochem.utils.sequence_utils as seq_utils


def get_bridging_oligos(target_melt_temp, sequences, plasmid_seq=None,
                        shuffle=False, reagent_concs=None):
    '''Designs dominos (bridging oligos) for LCR.'''
    num_sequences = len(sequences)
    orderings = sorted([list(a) for a in
                        set(itertools.permutations(range(num_sequences)))]) \
        if shuffle else [range(num_sequences)]

    if plasmid_seq is not None:
        sequences.append(plasmid_seq)
        for ordering in orderings:
            ordering.insert(0, num_sequences)
            ordering.append(num_sequences)

    pairs = []

    for ordering in orderings:
        pairs.extend(synbiochem.utils.pairwise(ordering))

    return [_get_bridge(sequences, num_sequences, pair, target_melt_temp,
                        reagent_concs)
            for pair in sorted(set(pairs))]


def _get_bridge(sequences, num_sequences, pair, target_melt_temp,
                reagent_concs=None):
    '''Get bridging oligo for pair of sequences.'''
    reverse, reverse_tm = _get_bridge_part(sequences[pair[0]], False,
                                           target_melt_temp, reagent_concs)
    forward, forward_tm = _get_bridge_part(sequences[pair[1]], True,
                                           target_melt_temp, reagent_concs)
    return ['SEQ' + str(pair[0] + 1) if pair[0] is not num_sequences else 'P',
            'SEQ' + str(pair[1] + 1) if pair[1] is not num_sequences else 'P',
            reverse, forward, str(reverse_tm), str(forward_tm),
            reverse + forward
            if reverse is not None and forward is not None else '']


def _get_bridge_part(sequence, forward, target_melt_temp, reagent_concs=None):
    '''Gets half of bridging oligo.'''
    for i in range(1, len(sequence)):
        subsequence = sequence[:(i + 1)] if forward else sequence[-(i + 1):]
        melting_temp = seq_utils.get_melting_temp(subsequence, None,
                                                  reagent_concs)

        if melting_temp > target_melt_temp:
            return subsequence, melting_temp

    return None, -float('inf')


def _read_input_file(filename):
    '''Reads input file.'''
    melting_temp = 70
    sequences = []
    plasmid_seq = None
    shuffle = False
    reagent_concs = {seq_utils.NA: 0.05,
                     seq_utils.K: 0,
                     seq_utils.TRIS: 0,
                     seq_utils.MG: 0.002,
                     seq_utils.DNTP: 1e-7}

    in_plasmid = False
    seq = ''

    with open(filename) as fle:
        for line in fle.read().splitlines():
            if line.startswith('MELTING_TEMP:'):
                melting_temp = float(line.replace('MELTING_TEMP:', '').strip())
            elif line.startswith('SHUFFLE:'):
                shuffle = line.replace('SHUFFLE:', '').strip() == 'True'
            elif line.startswith('NA:'):
                reagent_concs[seq_utils.NA] = \
                    float(line.replace('NA:', '').strip())
            elif line.startswith('K:'):
                reagent_concs[seq_utils.K] = \
                    float(line.replace('K:', '').strip())
            elif line.startswith('TRIS:'):
                reagent_concs[seq_utils.TRIS] = \
                    float(line.replace('TRIS:', '').strip())
            elif line.startswith('MG:'):
                reagent_concs[seq_utils.MG] = \
                    float(line.replace('MG:', '').strip())
            elif line.startswith('DNTP:'):
                reagent_concs[seq_utils.DNTP] = \
                    float(line.replace('DNTP:', '').strip())
            elif line.startswith('>'):
                if len(seq) > 0:
                    if in_plasmid:
                        plasmid_seq = seq
                    else:
                        sequences.append(seq)

                    seq = ''

                if line.startswith('>PLASMID'):
                    in_plasmid = True
            elif len(line) > 0 and not line.startswith('#') \
                    and seq is not None:
                seq += line.upper()

        if seq is not None:
            if in_plasmid:
                plasmid_seq = seq
                in_plasmid = False
            else:
                sequences.append(seq)

    return melting_temp, sequences, plasmid_seq, shuffle, reagent_concs


def main(argv):
    '''main method'''
    input_filename = argv[1]
    melting_temp, sequences, plasmid_seq, shuffle, reagent_concs = \
        _read_input_file(input_filename)
    results = get_bridging_oligos(melting_temp, sequences, plasmid_seq,
                                  shuffle, reagent_concs)

    result_file = open(os.path.splitext(input_filename)[0] +
                       '_result.xls', 'w+')

    for result in results:
        result_line = '\t'.join(result)
        print result_line
        result_file.write(result_line + '\n')

    result_file.close()

if __name__ == '__main__':
    main(sys.argv)

# get_bridging_oligos(70, ['saaaaae', 'sccccce', 'sggggge'], 'sttttte', True)
