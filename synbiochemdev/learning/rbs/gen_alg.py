'''
HolyGrail (c) University of Manchester 2015

HolyGrail is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments
import sys

from sbclearn import gen_alg

from synbiochemdev.learning.rbs import rbs_learn


def main(args):
    '''main method.'''

    # Get random peptides that match structure patterns from PDB:
    data, _ = rbs_learn.get_data(args[0])

    hyperparams = {
        # 'aa_props_filter': range(1, (2**holygrail.NUM_AA_PROPS)),
        # 'input_noise': [i / 10.0 for i in range(0, 10)],
        # 'hidden_noise': [i / 10.0 for i in range(0, 10)],
        'num_hidden_layers': range(1, 4),
        'num_nodes': range(10, 500, 10),
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

    predicter = rbs_learn.RBSRegression(data, float(args[1]), int(args[2]))

    learn_gen_alg = gen_alg.LearnGeneticAlgorithm(pop_size=int(args[3]),
                                                  predicter=predicter,
                                                  args=hyperparams,
                                                  retain=float(args[4]),
                                                  random_select=float(args[5]),
                                                  mutate=float(args[6]),
                                                  verbose=True)

    learn_gen_alg.run()


if __name__ == '__main__':
    main(sys.argv[1:])
