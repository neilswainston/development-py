'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from _collections import defaultdict

import random
import numpy

import matplotlib.pyplot as plt
import sbclearn.theanets.theanets_utils as theanets_utils


def learn(filename):
    '''Learn method.'''
    x_data, y_data, seqs = _get_data(filename)
    results = defaultdict(list)

    for _ in range(5):
        _learn(x_data, y_data, seqs, results)

    seqs = [val[0] for val in results.keys()]
    vals = [val[1] for val in results.keys()]
    preds = results.values()

    # The mean squared error:
    print('Mean squared error: %.2f'
          % numpy.mean([(x - y) ** 2
                        for x, y in zip(vals, [numpy.mean(pred)
                                               for pred in preds])]))

    for result in zip(seqs, vals,
                      [numpy.mean(pred) for pred in preds],
                      [numpy.std(pred) for pred in preds]):
        print '\t'.join([str(res) for res in result])

    _plot(vals, preds)


def _get_data(filename):
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

    return x_data[:len(y_data)], y_data, seqs


def _encode_x_data(x_data):
    '''Encodes x data.'''
    x_vals = {'A': (1, 0, 0, 0, 0),
              'C': (0, 1, 0, 0, 0),
              'G': (0, 0, 1, 0, 0),
              'T': (0, 0, 0, 1, 0),
              '-': (0, 0, 0, 0, 1)}

    return [val for nucl in x_data for val in x_vals[nucl]]


def _split_data(seqs, x_data, y_data, split=0.9):
    '''Split data.'''
    seqs_rand, x_data_rand, y_data_rand = \
        theanets_utils.randomise_order([seqs, x_data, y_data])

    # Split data into training and classifying:
    ind = int(split * len(x_data_rand))

    return seqs_rand[:ind], x_data_rand[:ind], \
        [[y] for y in y_data_rand[:ind]], \
        seqs_rand[ind:], x_data_rand[ind:], y_data_rand[ind:]


def _train(x_train, y_train, max_layers=3, max_nodes=100):
    '''Train neural network.'''
    regressor = theanets_utils.Regressor(x_train, y_train)
    hidden_layers = [random.randint(1, max_nodes)
                     for _ in range(random.randint(1, max_layers))]
    print str(hidden_layers)
    regressor.train(hidden_layers=hidden_layers)
    return regressor


def _learn(x_data, y_data, seqs, results):
    '''Learn method.'''
    _, x_train, y_train, seqs_val, x_val, y_val = \
        _split_data(seqs, x_data, y_data)

    regressor = _train(x_train, y_train)

    for tup in zip(*[seqs_val, y_val,
                     [val[0] for val in regressor.predict(x_val)]]):
        results[tup[:2]].append(tup[2])


def _plot(vals, preds):
    '''Plot results.'''
    plt.title('Prediction of limonene production from RBS seqs')
    plt.xlabel('Measured')
    plt.ylabel('Predicted')

    plt.errorbar(vals,
                 [numpy.mean(pred) for pred in preds],
                 yerr=[numpy.std(pred) for pred in preds],
                 fmt='o',
                 color='black')

    fit = numpy.poly1d(numpy.polyfit(vals,
                                     [numpy.mean(pred) for pred in preds], 1))
    plt.plot(vals, fit(vals), 'k')

    plt.xlim(0, 1.6)
    plt.ylim(0, 1.6)

    plt.show()


def main():
    '''main method.'''
    learn('rbs.txt')

if __name__ == '__main__':
    main()
