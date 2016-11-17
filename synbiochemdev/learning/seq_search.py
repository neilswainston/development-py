'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
# pylint: disable=too-few-public-methods
import operator
import random
import sys

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import matplotlib.pyplot
import numpy

from synbiochem.utils import sequence_utils
import sbclearn.theanets.theanets_utils as theanets_utils


class SequenceSearcher(object):
    '''Sequence activity search method.'''

    def __init__(self, seq_acts):
        self.__seq_acts = [list(val)
                           for val in sorted(seq_acts,
                                             key=operator.itemgetter(1))]
        self.__seq_acts.reverse()

    def get_seq_acts(self):
        '''Gets next generation of sequences to try.'''
        dissims = _calc_dissimilarities(self.__seq_acts[0][0],
                                        [val[0]
                                         for val in self.__seq_acts])

        self.__seq_acts = [vals + [dist] for vals, dist in zip(self.__seq_acts,
                                                               dissims)]
        self.__calc_pareto()

        return self.__seq_acts

    def __calc_pareto(self):
        '''Calculates the pareto front.'''
        # Start the Pareto frontier with the first value in the sorted list
        p_front = [self.__seq_acts[0][2]]
        self.__seq_acts[0].append(True)

        # Loop through the sorted list:
        for values in self.__seq_acts[1:]:
            pareto = False

            if values[2] >= p_front[-1]:  # Look for higher values of Y...
                # ...and add them to the Pareto frontier
                p_front.append(values[2])
                pareto = True

            values.append(pareto)


def _mutate(seq, max_mut_prob=0.1):
    '''Mutate sequence.'''
    mutant = ''
    mut_prob = random.random() * max_mut_prob

    for idx, _ in enumerate(seq):
        if random.random() < mut_prob:
            mutant += random.choice(['A', 'C', 'D', 'E', 'F', 'G', 'H',
                                     'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                                     'R', 'S', 'T', 'V', 'W', 'Y'])
        else:
            mutant += seq[idx]

    return mutant


def _get_random_seqs(length, num_seqs):
    '''Gets random sequences.'''
    seq = sequence_utils.get_random_aa(length)
    seqs = set([seq])

    while len(seqs) < num_seqs:
        mut = _mutate(seq)
        seqs.add(mut)

    return seq, seqs


def _calc_dissimilarities(ref_seq, seqs, matrix=None):
    '''Calculate sequence distances from reference sequence.'''
    if matrix is None:
        matrix = matlist.blosum62

    similarities = []
    gap_pen = -1000

    for seq in seqs:
        for alg in pairwise2.align.globalds(seq,
                                            ref_seq,
                                            matrix,
                                            gap_pen, gap_pen):

            similarities.append(alg[2])

    return [1 - dis for dis in _scale(similarities)]


def _scale(vals):
    '''Scales values.'''
    max_val = max(vals)
    min_val = min(vals)
    range_val = max_val - min_val
    return [(val - min_val) / range_val for val in vals]


def _plot_seq_acts(seq_acts):
    '''Plots sequence activity data.'''
    vals = zip(*seq_acts)
    pareto_vals = zip(*[p_vals for p_vals in seq_acts if p_vals[3]])

    matplotlib.pyplot.scatter(vals[2], vals[1],
                              c=['r' if val else 'b' for val in vals[3]])
    matplotlib.pyplot.plot(pareto_vals[2], pareto_vals[1], c='r')
    matplotlib.pyplot.xlabel('Sequence dis-similarity w.r.t. "best"')
    matplotlib.pyplot.ylabel('Activity')
    matplotlib.pyplot.show()


def _learn(sequences, activities):
    '''Attempt to learn sequence / activity relationship.'''
    # Convert sequences to inputs, based on amino acid properties:
    x_data = sequence_utils.get_aa_props(sequences)
    x_data, y_data = theanets_utils.randomise_order(x_data, activities)

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))

    y_train = [[y] for y in y_data[:ind]]
    regressor = theanets_utils.Regressor(x_data[:ind], y_train)

    regressor.train(hidden_layers=[1024])
    y_pred = regressor.predict(x_data[ind:])

    return regressor, y_data[ind:], y_pred


def _plot_activities(y_data, y_pred):
    '''Plots measured and predicted activities.'''
    for data, pred in zip(y_data, y_pred):
        print str(data) + '\t' + str(pred)

    matplotlib.pyplot.scatter(y_data, y_pred)
    matplotlib.pyplot.xlabel('Activity')
    matplotlib.pyplot.ylabel('Predicted activity')
    matplotlib.pyplot.show()

    print 'Mean delta: ' + str(numpy.mean([abs(j - i)
                                           for i, j in zip(y_data, y_pred)]))


def _output_pareto(seq_acts):
    '''Outputs pareto front members.'''
    pareto_seq_acts = [p_vals for p_vals in seq_acts if p_vals[3]]

    print '\t'.join(['Sequence', 'Activity',
                     'Sequence dis-similarity w.r.t. "best"'])

    for values in pareto_seq_acts:
        print '\t'.join([str(val) for val in values[0:3]])

    print '\n% pareto: {0:.2f}'.format(float(len(pareto_seq_acts)) /
                                       len(seq_acts) * 100.0)


def main(argv):
    '''main method.'''
    seq, seqs = _get_random_seqs(int(argv[0]), int(argv[1]))
    dissims = _calc_dissimilarities(seq, seqs)

    seq_acts = zip(seqs, [1 - dis + numpy.random.normal(scale=0.1)
                          for dis in dissims])

    seq_search = SequenceSearcher(seq_acts)
    seq_acts = seq_search.get_seq_acts()
    _plot_seq_acts(seq_acts)

    _, y_data, y_pred = _learn([val[0] for val in seq_acts],
                               [val[1] for val in seq_acts])

    _plot_activities(y_data, y_pred)
    _output_pareto(seq_acts)


if __name__ == '__main__':
    main(sys.argv[1:])
