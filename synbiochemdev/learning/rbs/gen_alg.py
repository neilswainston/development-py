'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
import sys

from synbiochemdev.learning.rbs import rbs_learn
import synbiochem.optimisation.gen_alg as gen_alg


class RbsGeneticAlgorithm(gen_alg.GeneticAlgorithm):
    '''Class to optimise parameters for Classifier using a GA.'''

    def __init__(self, pop_size, data, split, args,
                 retain=0.2, random_select=0.05, mutate=0.01, verbose=False):
        '''Constructor.'''
        super(RbsGeneticAlgorithm, self).__init__(pop_size, args, retain,
                                                  random_select, mutate,
                                                  verbose)

        self.__data = data
        self.__split = split

    def _fitness(self, individual):
        '''Determine the fitness of an individual.'''

        # Form hidden layers array:
        num_hidden_layers = individual.pop('num_hidden_layers', 1)
        activ_func = individual.pop('activ_func', 'relu')
        num_nodes = individual.pop('num_nodes', 128)
        individual['hidden_layers'] = [(num_nodes, activ_func)] * \
            num_hidden_layers

        cls = self.__runner.run(**individual)

        # Reform individual dict:
        individual.pop('hidden_layers')
        individual['num_hidden_layers'] = num_hidden_layers
        individual['activ_func'] = activ_func
        individual['num_nodes'] = num_nodes

        if self._verbose:
            print str(cls[3])
            print str(cls[4]) + '\t' + str(individual)
            print

        return 1 - cls[4]


def main(argv):
    '''main method.'''

    # Get random peptides that match structure patterns from PDB:
    data = rbs_learn.get_data('rbs.txt')

    args = {  # 'aa_props_filter': range(1, (2**holygrail.NUM_AA_PROPS)),
        # 'input_noise': [i / 10.0 for i in range(0, 10)],
        # 'hidden_noise': [i / 10.0 for i in range(0, 10)],
        # 'num_hidden_layers': range(1, 4),
        # 'num_nodes': range(100, 5000, 100),
        'activ_func': ['relu', 'prelu', 'lgrelu'],
        'learning_rate': [x / 10000.0 for x in range(1, 100)],
        'momentum': [x / 10.0 for x in range(0, 10)],
        'patience': range(1, 10),
        'min_improvement': [i / 1000.0 for i in range(1, 100)],
        # 'validate_every': range(1, 25),
        # 'batch_size': range(10, 50, 10),
        # 'hidden_dropout': [i * 0.1 for i in range(0, 10)],
        # 'input_dropout': [i * 0.1 for i in range(0, 10)]
    }

    rbs_gen_alg = RbsGeneticAlgorithm(pop_size=int(argv[1]),
                                      data=data,
                                      split=float(argv[5]),
                                      args=args,
                                      retain=float(argv[6]),
                                      random_select=float(argv[7]),
                                      mutate=float(argv[8]),
                                      verbose=True)

    rbs_gen_alg.run()


if __name__ == '__main__':
    main(sys.argv)
