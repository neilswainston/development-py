'''
Created on 13 Jun 2017

@author: neilswainston
'''
# pylint: disable=invalid-name
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main(args):
    '''main method.'''
    df = pd.read_table(args[0])
    df['WEIGHT'] = 1 / df['CONC'] ** 2
    x_range = range(-10, int(max(df['CONC'])))

    unweighted_fit = np.poly1d(np.polyfit(df['CONC'], df['PEAK AREA'], 1))
    weighted_fit = np.poly1d(np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                        w=df['WEIGHT']))

    plt.figure(1)
    plt.suptitle('Calibration curves')

    plt.subplot(121)
    _subplot(df, x_range, unweighted_fit, weighted_fit)

    plt.subplot(122)
    _subplot(df, x_range, unweighted_fit, weighted_fit)
    axes = plt.gca()
    axes.set_xlim([-1, 1])
    axes.set_ylim([-1, 1])
    plt.show()


def _subplot(df, x_range, unweighted_fit, weighted_fit):
    '''Make subplot.'''
    plt.xlabel(df.columns[0])
    plt.ylabel(df.columns[1])

    plt.scatter(df['CONC'], df['PEAK AREA'])

    unweighted, = plt.plot(x_range, unweighted_fit(x_range),
                           c='red', label='Unweighted: ' + str(unweighted_fit))
    weighted, = plt.plot(x_range, weighted_fit(x_range),
                         c='green', label='Weighted: ' + str(weighted_fit))

    plt.legend(handles=[unweighted, weighted])

    axes = plt.gca()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    axes.axvline(x=0, color='k')


if __name__ == '__main__':
    main(sys.argv[1:])
