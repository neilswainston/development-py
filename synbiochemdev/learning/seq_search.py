'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import matplotlib.pyplot
import operator
import random
import sys

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

from synbiochem.utils import sequence_utils


class SequenceSearcher(object):
    '''Sequence activity search method.'''

    def __init__(self, seq_act):
        self.__seq_act = [list(val)
                          for val in sorted(seq_act.items(),
                                            key=operator.itemgetter(1))]
        self.__seq_act.reverse()

    def get_next_gen(self):
        '''Gets next generation of sequences to try.'''
        distances = _calc_distances(self.__seq_act[0][0],
                                    [val[0] for val in self.__seq_act])

        self.__seq_act = [vals + [dist] for vals, dist in zip(self.__seq_act,
                                                              distances)]
        self.__calc_pareto()

        vals = zip(*self.__seq_act)
        pareto_vals = zip(*[p_vals for p_vals in self.__seq_act if p_vals[3]])

        matplotlib.pyplot.scatter(vals[2], vals[1],
                                  c=['r' if val else 'b' for val in vals[3]])
        matplotlib.pyplot.plot(pareto_vals[2], pareto_vals[1], c='r')
        matplotlib.pyplot.xlabel('Sequence dis-similarity w.r.t. "best"')
        matplotlib.pyplot.ylabel('Activity')
        matplotlib.pyplot.show()

        return [p_vals for p_vals in self.__seq_act if p_vals[3]]

    def __calc_pareto(self):
        '''Calculates the pareto front.'''
        # Start the Pareto frontier with the first value in the sorted list
        p_front = [self.__seq_act[0][2]]
        self.__seq_act[0].append(True)

        # Loop through the sorted list:
        for values in self.__seq_act[1:]:
            pareto = False

            if values[2] >= p_front[-1]:  # Look for higher values of Y...
                # ...and add them to the Pareto frontier
                p_front.append(values[2])
                pareto = True

            values.append(pareto)


def _mutate(seq, max_mut_prob=0.5):
    '''Mutate sequence.'''
    mutant = ''
    mut_prob = random.random() * max_mut_prob

    for idx in range(len(seq)):
        if random.random() < mut_prob:
            mutant += random.choice(['A', 'C', 'D', 'E', 'F', 'G', 'H',
                                     'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                                     'R', 'S', 'T', 'V', 'W', 'Y'])
        else:
            mutant += seq[idx]

    return mutant


def _calc_distances(ref_seq, seqs):
    '''Calculate sequence distances from best sequence.'''
    distances = []
    gap_pen = -1000
    best = None

    for seq in seqs:
        for alg in pairwise2.align.globalds(seq,
                                            ref_seq,
                                            matlist.blosum62,
                                            gap_pen, gap_pen):
            if best is None:
                best = alg[2]

            distances.append(1 - (alg[2] / best))

    return distances


def main(argv):
    '''main method.'''
    seq_act = {}
    seq = sequence_utils.get_random_aa(int(argv[0]))

    for _ in range(int(argv[1])):
        seq_act[_mutate(seq)] = random.random()

    seq_search = SequenceSearcher(seq_act)
    next_gen = seq_search.get_next_gen()

    print '\t'.join(['Sequence', 'Activity',
                     'Sequence dis-similarity w.r.t. "best"'])

    for values in next_gen:
        print '\t'.join([str(val) for val in values[0:3]])

if __name__ == '__main__':
    main(sys.argv[1:])
