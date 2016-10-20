'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math
import random
import sys

import numpy


def cent_limit_theorem(tests, lngth, min_val, max_val):
    '''Tests central limit theorem.'''
    means = []

    for _ in range(tests):
        rand_dist = get_rand_dist(lngth, min_val, max_val)
        mean_val = numpy.mean(rand_dist)
        print str(mean_val) + '\t' + str(numpy.std(rand_dist))
        means.append(mean_val)

    print '\t'.join([str(numpy.mean(means)),
                     str(numpy.std(means)),
                     str(numpy.std(means) * math.sqrt(lngth))])


def get_rand_dist(lngth, min_val, max_val):
    '''Returns random distribution.'''
    return [random.randrange(min_val, max_val) for _ in range(lngth)]


def main(argv):
    '''main method.'''
    return cent_limit_theorem(int(argv[0]), int(argv[1]), int(argv[2]),
                              int(argv[3]))

if __name__ == '__main__':
    main(sys.argv[1:])
