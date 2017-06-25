'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
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


def normal_cdf(x, mu=0, sigma=1):
    '''Returns normal cumulative distribution function.'''
    return (1 + math.erf((x - mu) / math.sqrt(2) / sigma)) / 2


def inverse_normal_cdf(p, mu=0.0, sigma=1.0, tolerance=0.00001):
    '''Find approx inverse using binary search.'''

    # If not standard, compute standard and rescale:
    if mu > 1e-4 or sigma != 1.0:
        return mu + sigma * inverse_normal_cdf(p, tolerance)

    lo_z = -10.0  # normal_cdf(-10) is ~0
    hi_z = 10.0  # normal_cdf(10) is ~1

    while hi_z - lo_z > tolerance:
        mid_z = (lo_z + hi_z) / 2
        mid_p = normal_cdf(mid_z)

        if mid_p < p:
            lo_z = mid_z
        elif mid_p > p:
            hi_z = mid_z
        else:
            break

    return mid_z


def main(argv):
    '''main method.'''
    return cent_limit_theorem(int(argv[0]), int(argv[1]), int(argv[2]),
                              int(argv[3]))


if __name__ == '__main__':
    main(sys.argv[1:])
