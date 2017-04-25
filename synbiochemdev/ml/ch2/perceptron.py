'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
import sys

from synbiochemdev.ml.ch2 import get_iris_data
import matplotlib.pyplot as plt
import numpy as np


class Perceptron(object):
    '''Perceptron classifier.'''

    def __init__(self, learn_rate=0.01, iters=10):
        '''Constructor.'''
        self.__learn_rate = learn_rate
        self.__iters = iters
        self.__model = {}

    def fit(self, X, y):
        '''Fit training data.'''
        self.__model['w'] = np.zeros(1 + X.shape[1])
        X = np.hstack((np.ones((X.shape[0], 1), dtype=X.dtype), X))
        errors = []

        for _ in range(self.__iters):
            curr_errs = 0

            for xi, target in zip(X, y):
                update = self.__learn_rate * (target - self.predict(xi))
                self.__model['w'] += update * xi
                curr_errs += int(update != 0.0)

            errors.append(curr_errs)

        return errors

    def predict(self, x):
        '''Return class label after unit step.'''
        return np.where(np.dot(x, self.__model['w']) >= 0.0, 1, -1)


def _plot_data(X):
    '''Plots data.'''
    plt.scatter(X[:50, 0], X[:50, 1], color='red', marker='o', label='setosa')
    plt.scatter(X[50:100, 0], X[50:100, 1], color='blue', marker='x',
                label='versicolor')

    plt.xlabel('sepal length')
    plt.ylabel('petal length')
    plt.legend(loc='upper left')

    plt.show()


def _plot_results(errors):
    '''Plots results.'''
    plt.plot(range(1, len(errors) + 1), errors, marker='o')
    plt.xlabel('Epochs')
    plt.ylabel('Number of misclassifications')
    plt.show()


def _run_ppn(learn_rate, iters, X, y):
    '''Runs perceptron.'''
    ppn = Perceptron(learn_rate=learn_rate, iters=iters)
    return ppn.fit(X, y)


def main(args):
    '''main method.'''
    data = get_iris_data()
    X = data.iloc[:, [0, 2]].values
    y = np.where(data.iloc[:, 4].values == 'Iris-setosa', -1, 1)

    # _plot_data(X)

    errors = _run_ppn(float(args[0]), int(args[1]), X, y)
    _plot_results(errors)

if __name__ == '__main__':
    main(sys.argv[1:])
