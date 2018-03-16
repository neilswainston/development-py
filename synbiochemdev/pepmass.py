'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
from collections import Counter, defaultdict
import random
import sys

from Bio.Seq import Seq
from synbiochem.utils import seq_utils

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyteomics import mass


def simulate(codon_templates, lngth, samples):
    '''Simulate random peptides.'''
    data = _get_data(codon_templates, lngth, samples)
    dct = {idx: pd.DataFrame(datum,
                             columns=['nucl_seq', 'prot_seq', 'mass'])
           for idx, datum in data.iteritems()}

    df = pd.concat(dct, axis=1, names=['barcode_num', 'fields'])
    df.to_csv('out.csv')

    df.dropna(inplace=True)
    print len(df)

    counter = Counter([tuple(sorted(vals))
                       for vals in df.loc[:, (slice(None), 'mass')].values])

    print len([val for val in counter.values() if val == 1])
    plot_hist(counter.values())


def _get_data(codon_templates, lngth, samples):
    '''Get data.'''
    data = defaultdict(list)

    for idx, codon_template in enumerate(codon_templates):
        for _ in range(samples):
            codons = [get_codon(codon_template) for _ in range(lngth)]
            nucl_seq = ''.join([nucl for codon in codons for nucl in codon])
            prot_seq = str(Seq(nucl_seq).translate())

            data[idx].append([nucl_seq,
                              prot_seq,
                              mass.calculate_mass(sequence=prot_seq)
                              if '*' not in prot_seq else float('NaN')])

    return data


def get_codon(codon_template):
    '''Get codon.'''
    return [random.choice(seq_utils.INV_NUCL_CODES[nucl])
            for nucl in codon_template]


def plot_hist(values):
    '''Plot histogram.'''
    plt.hist(values, bins=len(np.unique(values)))
    plt.xlabel('Number of monoisotopic peptides')
    plt.ylabel('Count')
    plt.title('Number of monoisotopic peptides histogram')
    plt.grid(True)
    plt.show()


def main(args):
    '''main method.'''
    simulate(args[2:], int(args[0]), int(args[1]))


if __name__ == '__main__':
    main(sys.argv[1:])
