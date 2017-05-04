'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
import sys
from synbiochemdev.ml.ch2 import get_iris_data, plot_epoch_progress, Classifier
import numpy as np


class Perceptron(Classifier):
    '''Perceptron classifier.'''

    def __init__(self, learn_rate=0.01, iters=10):
        '''Constructor.'''
        super(Perceptron, self).__init__(learn_rate, iters)

    def fit(self, X, y):
        '''Fit training data.'''
        self._model['w'] = np.zeros(1 + X.shape[1])
        X = np.hstack((np.ones((X.shape[0], 1), dtype=X.dtype), X))
        errors = []

        for _ in range(self._iters):
            curr_errs = 0

            for xi, target in zip(X, y):
                update = self._learn_rate * (target - self.predict(xi))
                self._model['w'] += update * xi
                curr_errs += int(update != 0.0)

            errors.append(curr_errs)

        return errors


def main(args):
    '''main method.'''
    data = get_iris_data()
    X = data.iloc[:, [0, 2]].values
    y = np.where(data.iloc[:, 4].values == 'Iris-setosa', -1, 1)

    # _plot_data(X)
    classifier = Perceptron(learn_rate=float(args[0]), iters=int(args[1]))
    epoch_progress = classifier.fit(X, y)
    plot_epoch_progress(epoch_progress, 'Errors')

if __name__ == '__main__':
    main(sys.argv[1:])
