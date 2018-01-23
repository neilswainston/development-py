'''
sbclearn (c) University of Manchester 2017

sbclearn is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
import sys

from sbclearn import gen_alg
from sbclearn.peptides import peptide_learn


def main(args):
    '''main method.'''

    # Get random peptides that match structure patterns from PDB:
    x_data, y_data = peptide_learn.get_data(args[0])

    hyperparams = {
        'input_noise': [i / 10.0 for i in range(0, 10)],
        'hidden_noise': [i / 10.0 for i in range(0, 10)],
        'num_hidden_layers': range(1, 4),
        'num_nodes': range(10, 500, 10),
        'activ_func': ['relu', 'prelu', 'lgrelu'],
        'learning_rate': [x / 10000.0 for x in range(1, 100)],
        'momentum': [x / 10.0 for x in range(0, 10)],
        'patience': range(1, 10),
        'min_improvement': [i / 1000.0 for i in range(1, 100)],
        'validate_every': range(1, 25),
        'batch_size': range(10, 50, 10),
        'hidden_dropout': [i * 0.1 for i in range(0, 10)],
        'input_dropout': [i * 0.1 for i in range(0, 10)]
    }

    learn_gen_alg = gen_alg.LearnGeneticAlgorithm(pop_size=int(args[3]),
                                                  x_data=x_data,
                                                  y_data=y_data,
                                                  test_size=float(args[1]),
                                                  tests=int(args[2]),
                                                  args=hyperparams,
                                                  retain=float(args[4]),
                                                  random_select=float(args[5]),
                                                  mutate=float(args[6]),
                                                  verbose=True)

    learn_gen_alg.run()


if __name__ == '__main__':
    main(sys.argv[1:])
