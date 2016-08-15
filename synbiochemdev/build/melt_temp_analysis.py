'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
from synbiochem.utils import sequence_utils

import matplotlib.pyplot as plt


def get_melt_temp(length=35, max_repeat_nuc=5, num=100):
    '''Gets melting temps.'''
    melt_tmps = []

    seqs = [sequence_utils.get_random_dna(length, max_repeat_nuc)
            for _ in range(num)]

    for idx, seq in enumerate(seqs):
        try:
            melt_tmps.extend([sequence_utils.get_melting_temp(seq, dna2,
                                                              strict=False)
                              for dna2 in seqs[idx + 1:]])
        except ZeroDivisionError, err:
            # Take no action
            print err

    return melt_tmps


def do_plot(data):
    '''Plots data.'''
    plt.hist(data, bins=100000)
    plt.xlabel('Tm')
    plt.ylabel('Count')
    plt.title('Tm histogram')
    plt.axis([-500, 100, 0, 1000])
    plt.grid(True)

    plt.savefig('tm.png')


def main():
    '''Main method.'''
    melt_tmps = get_melt_temp()
    do_plot(melt_tmps)

if __name__ == '__main__':
    main()
