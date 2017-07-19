'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import math

from sklearn.preprocessing.data import scale

import numpy as np


def de_mean_matrix(A):
    '''de-mean matrix.'''
    nr, nc = np.shape(A)
    column_means, _ = scale(A)
    return make_matrix(nr, nc, lambda i, j: A[i][j] - column_means[j])


def make_matrix(num_rows, num_cols, entry_fn):
    '''make matrix.'''
    return [[entry_fn(i, j) for j in range(num_cols)]
            for i in range(num_rows)]


def directional_variance(x, w):
    '''directional variance.'''
    return sum(directional_variance_i(x_i, w)
               for x_i in x)


def directional_variance_i(x_i, w):
    '''directional variance i.'''
    return dot(x_i, direction(w)) ** 2


def direction(w):
    '''direction.'''
    mag = magnitude(w)
    return [w_i / mag for w_i in w]


def magnitude(v):
    '''magnitude.'''
    return math.sqrt(sum_of_squares(v))


def sum_of_squares(v):
    '''sum of squares.'''
    return dot(v, v)


def dot(v, w):
    '''dot.'''
    return sum(v_i * w_i for v_i, w_i in zip(v, w))


def main():
    '''main method.'''


if __name__ == '__main__':
    main()
