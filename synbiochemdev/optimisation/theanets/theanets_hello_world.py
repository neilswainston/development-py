'''
Theanets "Hello, world!" - train a simple neural network for classifying
simple data and evaluate the results.

[Theanets](https://github.com/lmjohns3/theanets) allows to build and train
neural networks on top of the [Theano](https://github.com/Theano/Theano)
compiler.

The goal is to get familiar with theanets on some simple example. You can
modify this example bit by bit to work on more complex data and models.

In this example we generate some synthetic data (via scikit-learn) - two 2D
blobs with Gaussian distribution which are in addition linearly separable.
Thus any classification model should have no problem with such data.

We create a neural network with three layers - input, hidden and output - each
with two dimensions (2D featues, two classes). The input and hidden layer has
by default sigmoid activation, the output clasification layer has softmax
actiovation by default. The model is trained via the stochastic gradient
descent algorithm.

Finally the model is evaluated by functions provided by scikit-learn.
'''

# some utilities for command line interfaces
from collections import defaultdict
from functools import partial
from itertools import count

import climate
from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics import classification_report, confusion_matrix
import theanets

import numpy as np


# Deep neural networks on top of Theano:
def train(x_data, y_data, split=0.75, optimize='sgd', learning_rate=0.01,
          momentum=0.5):
    '''Train the network.'''

    # Check lengths of x_data and y_data are equal,
    # assume all tuples in x_data are of the same length.
    assert len(x_data) == len(y_data)

    ind = int(split * len(x_data))

    # Plain neural network with a single hidden layer:
    exp = theanets.Experiment(theanets.Classifier,
                              layers=(len(x_data[0]), 50, len(set(y_data))))

    # x_data = np.array(x_data, dtype=np.float32)
    y_map = _enumerate(y_data)
    y_data = np.array([y[1] for y in y_map], dtype=np.int32)

    exp.train((x_data[:ind], y_data[:ind]),
              (x_data[ind:], y_data[ind:]),
              optimize=optimize, learning_rate=learning_rate,
              momentum=momentum)

    return exp, dict(set(y_map))


def classify(exp, x_test, y_test, y_map):
    '''Classifies and analyses test data.'''
    y_test = np.array([y_map[y] for y in y_test], dtype=np.int32)
    y_pred = exp.network.classify(x_test)
    return y_map, y_pred, classification_report(y_test, y_pred), \
        confusion_matrix(y_test, y_pred)


def _enumerate(lst):
    '''Returns enumeration of supplied list.'''
    label_to_number = defaultdict(partial(next, count()))
    return [(item, label_to_number[item]) for item in lst]


def main():
    '''main method.'''
    climate.enable_default_logging()

    # -- generate some data --
    x_data, y_data = make_blobs(n_samples=1000, centers=5, n_features=3,
                                cluster_std=1.5, random_state=0)

    # convert the features and targets to the 32-bit format suitable for the
    # model
    x_data = x_data.astype(np.float32)
    # y_data = y_data.astype(np.int32)
    # x_data = x_data.tolist()
    y_data = [str(y) for y in y_data]

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))
    exp, y_map = train(x_data[:ind], y_data[:ind])

    for output in classify(exp, x_data[ind:], y_data[ind:], y_map):
        print output


if __name__ == '__main__':
    main()
