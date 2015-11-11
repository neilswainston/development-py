'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from itertools import chain
import random
import sys

from synbiochemdev.optimisation import theanets
from synbiochemdev.optimisation.theanets import randomise_order
import synbiochemdev.utils.structure_utils as struct_utils


# KD Hydrophobicity, EIIP, Helix, Sheet, Turn
AMINO_ACID_PROPS = {
    'A': [0.66, 0.336, 0.870, 0.399, 0.326],
    'R': [0.1, 0.707, 0.597, 0.518, 0.569],
    'N': [0.189, 0.123, 0.391, 0.276, 0.9],
    'D': [0.189, 0.9, 0.645, 0.276, 0.855],
    'C': [0.722, 0.625, 0.554, 0.608, 0.719],
    'Q': [0.189, 0.582, 0.691, 0.561, 0.59],
    'E': [0.189, 0.137, 0.9, 0.1, 0.403],
    'G': [0.464, 0.132, 0.213, 0.466, 0.9],
    'H': [0.216, 0.253, 0.691, 0.434, 0.569],
    'I': [0.9, 0.1, 0.589, 0.870, 0.1],
    'L': [0.838, 0.1, 0.814, 0.668, 0.251],
    'K': [0.153, 0.335, 0.755, 0.345, 0.61],
    'M': [0.669, 0.621, 0.828, 0.568, 0.263],
    'F': [0.749, 0.699, 0.683, 0.703, 0.263],
    'P': [0.358, 0.225, 0.1, 0.234, 0.883],
    'S': [0.429, 0.625, 0.323, 0.518, 0.842],
    'T': [0.438, 0.696, 0.391, 0.756, 0.576],
    'W': [0.42, 0.447, 0.622, 0.708, 0.576],
    'Y': [0.384, 0.427, 0.335, 0.746, 0.691],
    'V': [0.873, 0.136, 0.572, 0.9, 0.141]
}


def get_classif_data(sample_size, struct_patterns):
    '''Gets random data for classification analyses.'''
    return {struct_pattern: struct_utils.sample_seqs(sample_size,
                                                     struct_pattern)
            for struct_pattern in struct_patterns}


def get_regression_data(max_ids, num_samples, nmer_len):
    '''Gets random data for regression analyses.'''
    for pdb_id in struct_utils.get_pdb_ids(max_ids):
        try:
            for data in _sample_regression_data(pdb_id, num_samples,
                                                nmer_len):
                print '\t'.join([str(datum) for datum in data])
        except AssertionError, err:
            print err


def _sample_regression_data(pdb_id, num_samples, nmer_len):
    '''Samples learning data for deep learning.'''
    pdb_id, all_seqs, inpt, prox_output, phi_psi_output = \
        _get_regression_data(pdb_id)

    results = []

    num_aa_props = len(AMINO_ACID_PROPS['A'])

    for _ in range(num_samples):
        chn = int(random.random() * len(all_seqs))
        start = int(random.random() * (len(all_seqs[chn]) - nmer_len))
        rnge = range(start, start + nmer_len)
        result = [pdb_id,
                  ''.join([all_seqs[chn][i] for i in rnge]),
                  inpt[chn][start * num_aa_props:(start + nmer_len) *
                            num_aa_props],
                  [prox_output[chn][i]
                   for i in _get_sub_square_matrix(start, nmer_len,
                                                   len(all_seqs[chn]))],
                  phi_psi_output[chn][start * 2:(start + nmer_len) * 2]]

        results.append(result)

    return results


def _get_regression_data(pdb_id):
    '''Gets input/output for deep learning.'''
    all_sequences = struct_utils.get_sequences(pdb_id)
    all_proximities = struct_utils.calc_proximities(pdb_id)
    all_phi_psi_data = struct_utils.get_phi_psi_data(pdb_id)

    # Ensure data consistency in terms of data lengths:
    assert len(all_sequences) == len(all_proximities) == \
        len(all_phi_psi_data)
    assert all([len(all_sequences[idx]) == len(all_proximities[idx]) ==
                len(all_phi_psi_data[idx])
                for idx in range(len(all_sequences))])

    prox_output = [list(chain.from_iterable(proximities))
                   for proximities in all_proximities]

    phi_psi_output = [list(chain.from_iterable(lst))
                      for lst in all_phi_psi_data]

    return pdb_id, all_sequences, _get_input_data(all_sequences), \
        prox_output, phi_psi_output


def _get_sub_square_matrix(idx, lngth, size):
    '''Returns square submatrix of length lngth from square matrix of given
    size.'''
    return [((idx + r) * size) + idx + c for r in range(lngth)
            for c in range(lngth)]


def _get_input_data(all_sequences):
    '''Returns input data for machine-learning problems.'''
    return [list(chain.from_iterable([AMINO_ACID_PROPS[amino_acid]
                                      if amino_acid in AMINO_ACID_PROPS
                                      else [0.0, 0.0, 0.0]
                                      for amino_acid in sequences]))
            for sequences in all_sequences]


def main(argv):
    '''main method.'''
    classif_data = get_classif_data(int(argv[1]), argv[2:])

    x_data = _get_input_data([i for v in classif_data.values() for i in v])
    y_data = [i for k, v in classif_data.iteritems() for i in [k] * len(v)]

    x_data, y_data = randomise_order(x_data, y_data)

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))

    classifier = theanets.Classifier()
    classifier.train(x_data[:ind], y_data[:ind])

    for output in classifier.classify(x_data[ind:], y_data[ind:]):
        print output

if __name__ == '__main__':
    main(sys.argv)
