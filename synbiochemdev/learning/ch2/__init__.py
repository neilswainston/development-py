'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=unused-argument
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Classifier(object):
    '''Classifier.'''

    def __init__(self, learn_rate, iters):
        '''Constructor.'''
        self._learn_rate = learn_rate
        self._iters = iters
        self._model = {}

    def fit(self, X, y):
        '''Fit training data.'''
        return None

    def predict(self, x):
        '''Return class label after unit step.'''
        return np.where(np.dot(x, self._model['w']) >= 0.0, 1, -1)


def get_iris_data():
    '''Gets Iris dataset.'''
    return pd.read_csv('https://archive.ics.uci.edu/ml/'
                       'machine-learning-databases/iris/iris.data',
                       header=None)


def plot_iris_data(X):
    '''Plots data.'''
    plt.scatter(X[:50, 0], X[:50, 1], color='red', marker='o', label='setosa')
    plt.scatter(X[50:100, 0], X[50:100, 1], color='blue', marker='x',
                label='versicolor')

    plt.xlabel('sepal length')
    plt.ylabel('petal length')
    plt.legend(loc='upper left')

    plt.show()


def plot_epoch_progress(progress, ylabel):
    '''Plots progress by epoch.'''
    plt.plot(range(1, len(progress) + 1), progress, marker='o')
    plt.xlabel('Epochs')
    plt.ylabel(ylabel)
    plt.show()
