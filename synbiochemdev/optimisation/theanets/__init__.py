'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from collections import defaultdict
from functools import partial
from itertools import count

import climate
from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics import classification_report, confusion_matrix
import theanets

import numpy as np


class Classifier(object):
    '''Simple classifier in Theanets.'''

    def __init__(self, optimize='sgd', learning_rate=0.01,
                 momentum=0.5):
        self.__optimize = optimize
        self.__learning_rate = learning_rate
        self.__momentum = momentum
        self.__exp = None
        self.__y_map = None

    def train(self, x_data, y_data, split=0.75):
        '''Train the network.'''

        # Check lengths of x_data and y_data are equal,
        # assume all tuples in x_data are of the same length.
        assert len(x_data) == len(y_data)

        ind = int(split * len(x_data))

        # Plain neural network with a single hidden layer:
        self.__exp = theanets.Experiment(theanets.Classifier,
                                         layers=(len(x_data[0]), 50,
                                                 len(set(y_data))))

        x_data = np.array(x_data, dtype=np.float32)
        y_enum = _enumerate(y_data)
        y_data = np.array([y[1] for y in y_enum], dtype=np.int32)
        self.__y_map = dict(set(y_enum))

        self.__exp.train((x_data[:ind], y_data[:ind]),
                         (x_data[ind:], y_data[ind:]),
                         optimize=self.__optimize,
                         learning_rate=self.__learning_rate,
                         momentum=self.__momentum)

    def classify(self, x_test, y_test):
        '''Classifies and analyses test data.'''
        y_test = np.array([self.__y_map[y] for y in y_test], dtype=np.int32)
        y_pred = self.__exp.network.classify(x_test)
        inv_y_map = {v: k for k, v in self.__y_map.items()}
        return [inv_y_map[y] for y in y_pred], inv_y_map, \
            classification_report(y_test, y_pred), \
            confusion_matrix(y_test, y_pred)


def _enumerate(lst):
    '''Returns enumeration of supplied list.'''
    label_to_number = defaultdict(partial(next, count()))
    return [(item, label_to_number[item]) for item in lst]


def main():
    '''main method.'''
    climate.enable_default_logging()

    # Generate some data, convert to list of floats (inputs) and string 'names'
    # for output classifications:
    x_data, y_data = make_blobs(n_samples=1000, centers=5, n_features=3,
                                cluster_std=1.5, random_state=0)
    x_data = x_data.tolist()
    y_data = [str(unichr(y + ord('A'))) for y in y_data]

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))

    classifier = Classifier()
    classifier.train(x_data[:ind], y_data[:ind])

    for output in classifier.classify(x_data[ind:], y_data[ind:]):
        print output


if __name__ == '__main__':
    main()
