'''
Created on 5 Aug 2016

@author: neilswainston
'''
# pylint: disable=no-member
import matplotlib.pyplot
import numpy

from synbiochem.utils import sequence_utils
import holygrail.theanets_utils as theanets_utils


def _learn(sequences, activities):
    '''Attempt to learn sequence / activity relationship.'''
    # Convert sequences to inputs, based on amino acid properties:
    x_data = sequence_utils.get_aa_props(sequences)
    x_data, y_data = theanets_utils.randomise_order(x_data, activities)

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))

    y_train = [[y] for y in y_data[:ind]]
    regressor = theanets_utils.Regressor(x_data[:ind], y_train)

    regressor.train(hidden_layers=[1024])
    y_pred = regressor.predict(x_data[ind:])

    return regressor, y_data[ind:], y_pred


def main():
    '''main method.'''
    vals = []
    with open('Lynda_harris.txt', 'r') as fle:
        for line in fle.read().split('\r'):
            vals.append(line.split(','))

    sequences = [val[0] for val in vals]
    activities = [float(val[1]) for val in vals]

    _, y_data, y_pred = _learn(sequences, activities)

    for data, pred in zip(y_data, y_pred):
        print str(data) + '\t' + str(pred)

    matplotlib.pyplot.scatter(y_data, y_pred)
    matplotlib.pyplot.xlabel('Activity')
    matplotlib.pyplot.ylabel('Predicted activity')
    matplotlib.pyplot.show()

    print 'Mean delta: ' + str(numpy.mean([abs(j - i)
                                           for i, j in zip(y_data, y_pred)]))

if __name__ == '__main__':
    main()
