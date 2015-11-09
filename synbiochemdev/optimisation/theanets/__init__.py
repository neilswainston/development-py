'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from numpy.core import cumsum


def split(data, percentages):
    '''Splits data into percentage chunks.'''
    assert sum(percentages) == 1.0
    stops = [int(x * len(data)) for x in cumsum(percentages)]
    return [data[a:b] for a, b in zip([0]+stops, stops)]


print split(zip(range(20), range(20)), [0.2, 0.5, 0.3])
