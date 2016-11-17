'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import numpy

import matplotlib.pyplot as plt
import sbclearn.theanets.theanets_utils as theanets_utils


def _get_data():
    '''Gets data.'''
    x_data = []
    y_data = []

    with open('train.txt', 'rU') as infile:
        for line in infile:
            tokens = line.split('\t')
            x_data.append(tokens[0])
            y_data.append(float(tokens[1]))

    x_data = [vals for vals in zip(*[val
                                     for val in zip(*x_data)
                                     if len(set(val)) > 1])]

    x_data = [_encode_x_data(val) for val in x_data]

    return x_data, y_data


def _encode_x_data(x_data):
    '''Encodes x data.'''
    x_vals = {'A': (1, 0, 0, 0),
              'C': (0, 1, 0, 0),
              'G': (0, 0, 1, 0),
              'T': (0, 0, 0, 1)}

    return [val for nucl in x_data for val in x_vals[nucl]]


def _train(x_data, y_data):
    x_data_rand, y_data_rand = theanets_utils.randomise_order(x_data, y_data)

    # Split data into training and classifying:
    ind = int(0.9 * len(x_data_rand))

    y_train = [[y] for y in y_data_rand[:ind]]
    regressor = theanets_utils.Regressor(x_data_rand[:ind], y_train)

    regressor.train(hidden_layers=[50, 50])

    return regressor, x_data_rand, y_data_rand, ind


def main():
    '''main method.'''
    x_data, y_data = _get_data()

    all_vals = []
    all_pred = []

    for _ in range(100):
        regressor, x_data_rand, y_data_rand, ind = _train(x_data, y_data)

        all_pred.extend([val[0]
                         for val in regressor.predict(x_data_rand[ind:])])
        all_vals.extend(y_data_rand[ind:])

    plt.title('Prediction of limonene production from RBS seqs')
    plt.xlabel('Measured')
    plt.ylabel('Predicted')

    plt.scatter(all_vals, all_pred)

    fit = numpy.poly1d(numpy.polyfit(all_vals, all_pred, 1))
    plt.plot(all_vals, fit(all_vals), 'r')

    plt.xlim(0, 1.6)
    plt.ylim(0, 1.6)

    plt.show()


if __name__ == '__main__':
    main()
