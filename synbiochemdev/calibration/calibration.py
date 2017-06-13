'''
Created on 13 Jun 2017

@author: neilswainston
'''
# pylint: disable=invalid-name
import sys

from scipy import stats

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main(args):
    '''main method.'''
    df = pd.read_table(args[0])
    x_range = range(-10, int(max(df['CONC'])))

    weightings = {
        'Unweighted': np.polyfit(df['CONC'], df['PEAK AREA'], 1),
        'Weighted (1/x^0.5)': np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                         w=(1 / df['CONC'] ** 0.5)),
        'Weighted (1/x)': np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                     w=(1 / df['CONC'])),
        'Weighted (1/x^2)': np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                       w=(1 / df['CONC'] ** 2)),
    }

    plt.figure(1)
    plt.suptitle('Calibration curves')

    plt.subplot(121)
    _subplot(df, x_range, weightings)

    plt.subplot(122)
    _subplot(df, x_range, weightings)
    axes = plt.gca()
    axes.set_xlim([-1, 12])
    axes.set_ylim([-0.025, 0.25])
    plt.show()


def _subplot(df, x_range, weightings):
    '''Make subplot.'''
    plt.xlabel(df.columns[0])
    plt.ylabel(df.columns[1])

    plt.scatter(df['CONC'], df['PEAK AREA'])

    handles = []

    for label, weight in weightings.iteritems():
        weight = np.poly1d(weight)
        _, _, r_value, _, _ = \
            stats.linregress(df['PEAK AREA'], weight(df['CONC']))
        # label += ': R2=%.3f' % r_value
        ret = plt.plot(x_range, weight(x_range), label=label)

        handles.append(ret[0])

    plt.legend(handles=handles)

    axes = plt.gca()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    axes.axvline(x=0, color='k')


if __name__ == '__main__':
    main(sys.argv[1:])
