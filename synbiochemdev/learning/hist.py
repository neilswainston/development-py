'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
import uuid

import matplotlib.pyplot as plt
import numpy as np


def do_plot(data):
    '''Plots data.'''
    for datum, colour in data:
        plt.hist(datum, 50, facecolor=colour, alpha=0.75)

    plt.xlabel('Activity')
    plt.ylabel('Count')
    plt.title('Activity histogram')
    plt.axis([0, 8, 0, 10000])
    plt.grid(True)

    plt.savefig(str(uuid.uuid4()) + '.png')


def main():
    '''main method.'''
    poor = (1 + 0.2 * np.random.randn(100000), 'red')
    good = (2 + 0.4 * np.random.randn(100000), 'orange')
    great = (5 + 1.0 * np.random.randn(100000), 'green')

    do_plot([poor])
    do_plot([poor, good])
    do_plot([poor, good, great])


if __name__ == '__main__':
    main()
