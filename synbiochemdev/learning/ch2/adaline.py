'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
import sys

from synbiochemdev.learning.ch2 import get_iris_data, plot_epoch_progress, \
    Classifier
import numpy as np


class Adaline(Classifier):
    '''Adaline classifier.'''

    def __init__(self, learn_rate=0.01, iters=10):
        '''Constructor.'''
        super(Adaline, self).__init__(learn_rate, iters)

    def fit(self, X, y):
        '''Fit training data.'''
        self._model['w'] = np.zeros(1 + X.shape[1])
        X = np.hstack((np.ones((X.shape[0], 1), dtype=X.dtype), X))
        cost = []

        for _ in range(self._iters):
            output = np.dot(X, self._model['w'])
            errors = (y - output)
            self._model['w'] += self._learn_rate * X.T.dot(errors)
            cost.append((errors**2).sum() / 2.0)

        return cost


def main(args):
    '''main method.'''
    data = get_iris_data()
    X = data.iloc[:, [0, 2]].values
    y = np.where(data.iloc[:, 4].values == 'Iris-setosa', -1, 1)

    # _plot_data(X)
    classifier = Adaline(learn_rate=float(args[0]), iters=int(args[1]))
    epoch_progress = classifier.fit(X, y)
    plot_epoch_progress(epoch_progress, 'Cost')


if __name__ == '__main__':
    main(sys.argv[1:])
