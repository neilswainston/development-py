'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
import sys

from synbiochem.utils import chem_utils

from sbclearn.theanets.theanets_utils import Regressor
import pandas as pd
import sbclearn


def main(args):
    '''main method.'''
    data = pd.read_csv(args[0])
    split = float(args[1])

    X = [chem_utils.get_fingerprint(smiles, radius=8).tolist()
         for smiles in data.SMILES.tolist()]
    y = [[val] for val in data.iloc[:, 6].values]

    hyperparams = {
        # 'input_noise': [i / 10.0 for i in range(0, 10)],
        # 'hidden_noise': [i / 10.0 for i in range(0, 10)],
        'activ_func': 'relu',
        'learning_rate': 0.001,
        'momentum': 0.5,
        'patience': 3,
        'min_improvement': 0.1,
        # 'validate_every': range(1, 25),
        # 'batch_size': 10,
        # 'hidden_dropout': [i * 0.1 for i in range(0, 10)],
        # 'input_dropout': [i * 0.1 for i in range(0, 10)]
    }

    for _ in range(25):
        X, y = sbclearn.randomise_order((X, y))
        regressor = Regressor(X[:int(len(X) * split)],
                              y[:int(len(X) * split)])

        regressor.train(hidden_layers=[250, 100], hyperparams=hyperparams)
        y_preds, _ = regressor.predict(X[int(len(X) * split):])

        for val, pred in zip(y[int(len(X) * split):], y_preds):
            print '\t'.join([str(val[0]), str(pred)])

if __name__ == '__main__':
    main(sys.argv[1:])
