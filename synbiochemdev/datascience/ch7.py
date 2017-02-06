'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import math
import sys

from synbiochemdev.datascience.ch6 import inverse_normal_cdf, normal_cdf


# Normal cdf is the prob that teh variable is below a threshold:
normal_prob_below = normal_cdf


def normal_prob_above(lo, mu=0, sigma=1):
    '''It's above the threshold if not below...'''
    return 1 - normal_prob_below(lo, mu, sigma)


def normal_prob_between(lo, hi, mu=0, sigma=1):
    '''It's between if less than high, but not less than low...'''
    return normal_prob_below(hi, mu, sigma) - normal_prob_below(lo, mu, sigma)


def normal_prob_outside(lo, hi, mu=0, sigma=1):
    '''It's outside if not between...'''
    return 1 - normal_prob_between(lo, hi, mu, sigma)


def normal_upper_bound(prob, mu=0, sigma=1):
    '''Returns the z for which P(Z <= z) = probability...'''
    return inverse_normal_cdf(prob, mu, sigma)


def normal_lower_bound(prob, mu=0, sigma=1):
    '''Returns the x for which P(Z >= z) = probability...'''
    return inverse_normal_cdf(1 - prob, mu, sigma)


def normal_two_sided_bounds(prob, mu=0, sigma=1):
    '''Returns the symmetric (about the mean) bounds that contains the
    specified probability.'''
    tail_prob = (1 - prob) / 2

    # Upper bound should have tail_prob above it:
    upper_bound = normal_lower_bound(tail_prob, mu, sigma)

    # Lower bound should have tail_prob below it:
    lower_bound = normal_upper_bound(tail_prob, mu, sigma)

    return lower_bound, upper_bound


def normal_approximation_to_binomial(n, p):
    '''Finds mu and sigma corresponding to Binomial(n,p).'''
    mu = p * n
    sigma = math.sqrt(p * (1 - p) * n)
    return mu, sigma


def main(args):
    '''main method.'''
    mu_0, sigma_0 = normal_approximation_to_binomial(int(args[0]),
                                                     float(args[1]))

    print normal_two_sided_bounds(0.95, mu_0, sigma_0)

if __name__ == '__main__':
    main(sys.argv[1:])
