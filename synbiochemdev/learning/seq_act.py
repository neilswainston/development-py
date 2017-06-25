'''
Created on 5 Aug 2016

@author: neilswainston
'''
# pylint: disable=no-member
from synbiochem.utils import seq_utils
import matplotlib.pyplot
import numpy

from sbclearn import theanets
import sbclearn


def _learn(sequences, activities):
    '''Attempt to learn sequence / activity relationship.'''
    # Convert sequences to inputs, based on amino acid properties:
    x_data = sbclearn.get_aa_props(sequences)
    x_data, y_data = sbclearn.randomise_order((x_data, activities))

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))

    y_train = [[y] for y in y_data[:ind]]
    regressor = theanets.theanets_utils.Regressor(x_data[:ind], y_train)

    regressor.train(hidden_layers=[1024, 1024])
    y_pred = regressor.predict(x_data[ind:])

    return regressor, y_data[ind:], y_pred


def _plot(y_data, y_pred):
    '''Plots data.'''
    matplotlib.pyplot.scatter(y_data, y_pred)
    matplotlib.pyplot.xlabel('Activity')
    matplotlib.pyplot.ylabel('Predicted activity')
    matplotlib.pyplot.show()

    for data, pred in zip(y_data, y_pred):
        print str(data) + '\t' + str(pred)

    print 'Mean delta: ' + str(numpy.mean([abs(j - i)
                                           for i, j in zip(y_data, y_pred)]))


def main():
    '''main method.'''
    mao = seq_utils.get_uniprot_values(['P46882'], ['sequence'])
    mao_seq = mao['P46882']['Sequence']

    with open('raw.csv', 'r') as raw, open('seq_act.txt', 'w') as seq_act:
        for line in raw.read().split('\r'):
            tokens = line.split(',')
            offset = int(tokens[0]) - 1
            seq_act.write(mao_seq[:offset] +
                          tokens[1] +
                          mao_seq[offset + len(tokens[1]):] +
                          ',' + tokens[2] + '\r')

    vals = []
    with open('seq_act.txt', 'r') as fle:
        for line in fle.read().strip().split('\r'):
            vals.append(line.split(','))

    sequences = [val[0] for val in vals]
    activities = [float(val[1]) for val in vals]

    _, y_data, y_pred = _learn(sequences, activities)
    _plot(y_data, y_pred)


if __name__ == '__main__':
    main()
