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


def fit_curves(df):
    '''Perform analysis.'''
    curves = []

    for coeff in np.arange(0, 3, 1):
        curve = np.poly1d(np.polyfit(df['CONC'], df['PEAK AREA'], 1,
                                     w=(1 / df['CONC'] ** coeff)))
        ae = np.abs((df['PEAK AREA'] - curve(df['CONC'])) / df['PEAK AREA'])
        curves.append([coeff, curve, np.mean(ae)])

    return curves


def map_peak_areas(curves, peak_areas):
    '''Map hypothetical peak area data.'''
    concs = {}

    for peak_area in peak_areas:

        concs[peak_area] = [(curve[1] - peak_area).roots[0]
                            for curve in curves]

        print '%f\t%s' % (peak_area, str(concs[peak_area]))

    return concs


def plot(idx, max_idx, df, curves):
    '''Plot data.'''
    x_plot_range = sorted(set(df['CONC'].tolist()))[1] * 1.2

    plt.subplot(max_idx, 2, idx * 2 + 1)
    _subplot(df, curves)

    plt.subplot(max_idx, 2, idx * 2 + 2)
    _subplot(df, curves)

    axes = plt.gca()
    axes.set_xlim([-x_plot_range, x_plot_range])
    axes.set_ylim([-1000, 1000])


def _subplot(df, curves):
    '''Make subplot.'''
    x_plot_range = sorted(set(df['CONC'].tolist()))[1] * 1.2
    x_range = range(-int(math.ceil(x_plot_range)), int(max(df['CONC'])))

    plt.xlabel(df.columns[0])
    plt.ylabel(df.columns[1])

    plt.scatter(df['CONC'], df['PEAK AREA'], marker='x', color='k')

    handles = []

    for curve in curves:
        label = '(1/x^%.1f): mae=%.3f, curve: %s' % (
            curve[0], curve[2], curve[1])
        ret = plt.plot(x_range, curve[1](x_range), label=label, linestyle=':')
        handles.append(ret[0])

    plt.legend(handles=handles)

    axes = plt.gca()
    axes.grid(True, which='both')
    axes.axhline(y=0, color='k')
    axes.axvline(x=0, color='k')


def main(args):
    '''main method.'''
    plt.figure(1)
    plt.suptitle('Calibration curves')

    for idx, filename in enumerate(args):
        df = pd.read_table(filename)
        curves = fit_curves(df)
        plot(idx, len(args), df, curves)

    plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
