'''
Created on 13 Jun 2017

@author: neilswainston
'''
# pylint: disable=invalid-name
import math
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def analyse(idx, max_idx, filename):
    '''Perform analysis.'''
    df = pd.read_table(filename)
    x_plot_range = sorted(set(df['CONC'].tolist()))[1] * 1.2
    x_range = range(-int(math.ceil(x_plot_range)), int(max(df['CONC'])))

    weightings = [
        ['Unweight', np.polyfit(df['CONC'], df['PEAK AREA'], 1)],
        ['Weight (1/x^0.5)', np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                        w=(1 / df['CONC'] ** 0.5))],
        ['Weight (1/x)', np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                    w=(1 / df['CONC']))],
        ['Weight (1/x^2)', np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                      w=(1 / df['CONC'] ** 2))]
    ]

    plt.subplot(max_idx, 2, idx * 2 + 1)
    plt.title = filename
    _subplot(df, x_range, weightings)

    plt.subplot(max_idx, 2, idx * 2 + 2)
    plt.title = filename
    _subplot(df, x_range, weightings)

    axes = plt.gca()
    axes.set_xlim([-x_plot_range, x_plot_range])
    axes.set_ylim([-0.025, 0.25])


def main(args):
    '''main method.'''
    plt.figure(1)
    plt.suptitle('Calibration curves')

    for idx, filename in enumerate(args):
        analyse(idx, len(args), filename)

    plt.show()


def _subplot(df, x_range, weightings):
    '''Make subplot.'''
    plt.xlabel(df.columns[0])
    plt.ylabel(df.columns[1])

    plt.scatter(df['CONC'], df['PEAK AREA'])

    handles = []

    for weight in weightings:
        wgt = np.poly1d(weight[1])
        ae = np.abs((df['PEAK AREA'] - wgt(df['CONC'])) / df['PEAK AREA'])
        label = weight[0] + ': mae=%.3f' % np.mean(ae)
        ret = plt.plot(x_range, wgt(x_range), label=label)

        handles.append(ret[0])

    plt.legend(handles=handles)

    axes = plt.gca()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    axes.axvline(x=0, color='k')


if __name__ == '__main__':
    main(sys.argv[1:])
