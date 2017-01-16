'''
Created on 5 Aug 2016

@author: neilswainston
'''
# pylint: disable=no-member
import matplotlib.pyplot
import numpy

import holygrail.theanets_utils as theanets_utils


def _learn(inputs, outputs):
    '''Attempt to learn sequence / activity relationship.'''
    x_data, y_data = theanets_utils.randomise_order(inputs, outputs)

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))

    y_train = [[y] for y in y_data[:ind]]
    regressor = theanets_utils.Regressor(x_data[:ind], y_train)

    regressor.train(hidden_layers=[1])
    y_pred = regressor.predict(x_data[ind:])

    return regressor, y_data[ind:], y_pred


def main():
    '''main method.'''
    vals = []
    with open('houses.csv', 'r') as fle:
        for line in fle.read().split('\r'):
            vals.append(line.split(','))

    inputs = [[float(val[0]), float(val[1])] for val in vals]
    outputs = [float(val[2]) for val in vals]

    _, y_data, y_pred = _learn(inputs, outputs)

    for data, pred in zip(y_data, y_pred):
        print str(data) + '\t' + str(pred)

    matplotlib.pyplot.scatter(y_data, y_pred)
    matplotlib.pyplot.xlabel('Price')
    matplotlib.pyplot.ylabel('Predicted price')
    matplotlib.pyplot.show()

    print 'Mean delta: ' + str(numpy.mean([abs(j - i)
                                           for i, j in zip(y_data, y_pred)]))

if __name__ == '__main__':
    main()
