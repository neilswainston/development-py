'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-locals
import math

from synbiochemdev.datascience.ch6 import inverse_normal_cdf, normal_cdf
import matplotlib.pyplot as plt


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


def two_sided_p_value(x, mu, sigma):
    '''Gets two-sided p-value.'''
    if x >= mu:
        # Tail is greater than x:
        return 2 * normal_prob_above(x, mu, sigma)

    # Tail is less than x:
    return 2 * normal_prob_below(x, mu, sigma)


def main():
    '''main method.'''
    from numpy.random import random

    # Distribution if n=1000, p=0.5:
    mu_0, sigma_0 = normal_approximation_to_binomial(1000, 0.5)
    print mu_0, sigma_0

    # 95% bounds:
    lo_0, hi_0 = normal_two_sided_bounds(0.95, mu_0, sigma_0)
    print lo_0, hi_0

    # Distribution if n=1000, p=0.55:
    mu_1, sigma_1 = normal_approximation_to_binomial(1000, 0.55)
    print mu_1, sigma_1

    x_data = range(400, 600)
    y_data_0 = [normal_cdf(x + 0.5, mu_0, sigma_0) -
                normal_cdf(x - 0.5, mu_0, sigma_0)
                for x in x_data]
    y_data_1 = [normal_cdf(x + 0.5, mu_1, sigma_1) -
                normal_cdf(x - 0.5, mu_1, sigma_1)
                for x in x_data]

    # Type 2 (false positive) probability if x is in original interval:
    # Two-sided-test:
    type_2_prob_0 = normal_prob_between(lo_0, hi_0, mu_0, sigma_0)
    power_0 = 1 - type_2_prob_0
    print power_0

    type_2_prob_1 = normal_prob_between(lo_0, hi_0, mu_1, sigma_1)
    power_1 = 1 - type_2_prob_1
    print power_1

    # One-sided test:
    hi_2 = normal_upper_bound(0.95, mu_0, sigma_0)
    print hi_2

    type_2_prob_2 = normal_prob_below(hi_2, mu_1, sigma_1)
    power_2 = 1 - type_2_prob_2
    print power_2

    # p-value:
    print two_sided_p_value(529.5, mu_0, sigma_0)

    # p-value validation through simulation:
    extreme_val_count = 0

    for _ in range(100000):
        num_heads = sum(1 if random() < 0.5 else 0 for _ in range(1000))

        if num_heads >= 530 or num_heads <= 470:
            extreme_val_count += 1

    print extreme_val_count / 100000.0

    upper_p_val = normal_prob_above
    # lower_p_val = normal_prob_below

    print upper_p_val(524.5, mu_0, sigma_0)
    print upper_p_val(526.5, mu_0, sigma_0)

    # Plots:
    plt.plot(x_data, y_data_0)
    plt.plot(x_data, y_data_1)
    plt.xlabel('x')
    plt.ylabel('y price')
    plt.show()


if __name__ == '__main__':
    main()
