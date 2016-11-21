'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from _collections import defaultdict

import numpy

import matplotlib.pyplot as plt
import sbclearn.theanets.theanets_utils as theanets_utils


def _get_data(filename):
    '''Gets data.'''
    x_data = []
    y_data = []
    test_seqs = []

    with open(filename, 'rU') as infile:
        for line in infile:
            tokens = line.strip().split('\t')
            x_data.append(tokens[0])

            if len(tokens) > 1:
                y_data.append(float(tokens[1]))
            else:
                test_seqs.append(tokens[0])

    x_data = [vals for vals in zip(*[val
                                     for val in zip(*x_data)
                                     if len(set(val)) > 1])]

    x_data = [_encode_x_data(val) for val in x_data]

    return x_data[:len(y_data)], y_data, x_data[len(y_data):], test_seqs


def _encode_x_data(x_data):
    '''Encodes x data.'''
    x_vals = {'A': (1, 0, 0, 0, 0),
              'C': (0, 1, 0, 0, 0),
              'G': (0, 0, 1, 0, 0),
              'T': (0, 0, 0, 1, 0),
              '-': (0, 0, 0, 0, 1)}

    return [val for nucl in x_data for val in x_vals[nucl]]


def _split_data(x_data, y_data, split=0.9):
    '''Split data.'''
    x_data_rand, y_data_rand = theanets_utils.randomise_order(x_data, y_data)

    # Split data into training and classifying:
    ind = int(split * len(x_data_rand))

    x_train = x_data_rand[:ind]
    y_train = [[y] for y in y_data_rand[:ind]]
    x_val = x_data_rand[ind:]
    y_val = y_data_rand[ind:]

    return x_train, y_train, x_val, y_val


def _train(x_train, y_train):
    '''Train neural network.'''
    regressor = theanets_utils.Regressor(x_train, y_train)
    regressor.train(hidden_layers=[50, 25, 10])
    return regressor


def main():
    '''main method.'''
    x_data, y_data, x_test, test_seqs = _get_data('rbs.txt')
    val_vals = []
    val_pred = []
    test_pred = []

    for _ in range(50):
        x_train, y_train, x_val, y_val = _split_data(x_data, y_data)
        regressor = _train(x_train, y_train)

        val_pred.extend([val[0]
                         for val in regressor.predict(x_val)])
        val_vals.extend(y_val)

        test_pred.append([val[0]
                          for val in regressor.predict(x_test)])

    plt.title('Prediction of limonene production from RBS seqs')
    plt.xlabel('Measured')
    plt.ylabel('Predicted')

    results = defaultdict(list)

    for val, pred in zip(val_vals, val_pred):
        results[val].append(pred)

    plt.errorbar(results.keys(),
                 [numpy.mean(vals) for vals in results.values()],
                 yerr=[numpy.std(vals) for vals in results.values()],
                 fmt='o',
                 color='black')

    test_pred = zip(*test_pred)

    plt.errorbar([numpy.mean(vals) for vals in test_pred],
                 [numpy.mean(vals) for vals in test_pred],
                 yerr=[numpy.std(vals) for vals in test_pred],
                 fmt='o',
                 color='red')

    fit = numpy.poly1d(numpy.polyfit(val_vals, val_pred, 1))
    plt.plot(val_vals, fit(val_vals), 'k')

    plt.xlim(0, 1.6)
    plt.ylim(0, 1.6)

    # The mean squared error:
    print('Mean squared error: %.2f'
          % numpy.mean([(x - y) ** 2 for x, y in zip(val_vals, val_pred)]))

    for seq, pred, s_d in zip(test_seqs,
                              [numpy.mean(vals) for vals in test_pred],
                              [numpy.std(vals) for vals in test_pred]):
        print seq + '\t' + str(pred) + '\t' + str(s_d)

    plt.show()


if __name__ == '__main__':
    main()
