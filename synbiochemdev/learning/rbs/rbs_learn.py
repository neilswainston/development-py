'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
from _collections import defaultdict

import numpy

import matplotlib.pyplot as plt
import sbclearn.theanets.theanets_utils as theanets_utils


class RBSRegression(object):
    '''Class to perform prediction of production from RBS sequence.'''

    def __init__(self, data, split=0.9, tests=1):
        self.__data = data
        self.__split = split
        self.__tests = tests

    def predict(self, hidden_layers=None, hyperparams=None):
        '''Learn method.'''
        results = defaultdict(list)

        for _ in range(self.__tests):
            _, error = self.__learn(hidden_layers, hyperparams, results)

        return error, results

    def __learn(self, hidden_layers, hyperparams, results):
        '''Learn method.'''
        x_train, y_train, x_val, y_val = \
            theanets_utils.split_data(self.__data, self.__split)

        regressor = theanets_utils.Regressor(x_train, y_train)
        regressor.train(hidden_layers=hidden_layers, hyperparams=hyperparams)

        return regressor.predict(x_val, y_val, results=results)


def get_data(filename):
    '''Gets data.'''
    x_data = []
    y_data = []

    with open(filename, 'rU') as infile:
        for line in infile:
            tokens = line.strip().split('\t')
            x_data.append(tokens[0])

            if len(tokens) > 1:
                y_data.append(float(tokens[1]))

    x_data = [''.join(vals) for vals in zip(*[val
                                              for val in zip(*x_data)
                                              if len(set(val)) > 1])]
    seqs = x_data
    x_data = [_encode_x_data(val) for val in x_data]

    return [x_data, y_data], seqs


def _encode_x_data(x_data):
    '''Encodes x data.'''
    x_vals = {'A': (1, 0, 0, 0, 0),
              'C': (0, 1, 0, 0, 0),
              'G': (0, 0, 1, 0, 0),
              'T': (0, 0, 0, 1, 0),
              '-': (0, 0, 0, 0, 1)}

    return [val for nucl in x_data for val in x_vals[nucl]]


def _output(error, results):
    '''Output results.'''
    print 'Mean squared error: %.3f' % error

    for result in zip(results.keys(),
                      [numpy.mean(pred) for pred in results.values()],
                      [numpy.std(pred) for pred in results.values()]):
        print '\t'.join([str(res) for res in result])

    _plot(results)


def _plot(results):
    '''Plot results.'''
    plt.title('Prediction of limonene production from RBS seqs')
    plt.xlabel('Measured')
    plt.ylabel('Predicted')

    plt.errorbar(results.keys(),
                 [numpy.mean(pred) for pred in results.values()],
                 yerr=[numpy.std(pred) for pred in results.values()],
                 fmt='o',
                 color='black')

    fit = numpy.poly1d(numpy.polyfit(results.keys(),
                                     [numpy.mean(pred)
                                      for pred in results.values()], 1))
    plt.plot(results.keys(), fit(results.keys()), 'k')

    plt.xlim(0, 1.6)
    plt.ylim(0, 1.6)

    plt.show()


def main():
    '''main method.'''
    data, _ = get_data('rbs.txt')
    rbs_reg = RBSRegression(data)
    error, results = rbs_reg.predict([10, 10, 10])
    _output(error, results)

if __name__ == '__main__':
    main()
