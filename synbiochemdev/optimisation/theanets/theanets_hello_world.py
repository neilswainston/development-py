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
import climate
# deep neural networks on top of Theano
import theanets
import numpy as np
from sklearn.datasets.samples_generator import make_blobs
from sklearn.metrics import classification_report, confusion_matrix


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

    exp.train((x_data[:ind], y_data[:ind]),
              (x_data[ind:], y_data[ind:]),
              optimize=optimize, learning_rate=learning_rate,
              momentum=momentum)

    return exp


def classify(exp, x_test, y_test):
    '''Classifies and analyses test data.'''
    y_pred = exp.network.classify(x_test)
    return y_pred, classification_report(y_test, y_pred), \
        confusion_matrix(y_test, y_pred)


def main():
    '''main method.'''
    climate.enable_default_logging()

    # -- generate some data --
    x_data, y_data = make_blobs(n_samples=1000, centers=5, n_features=3,
                                cluster_std=1.5, random_state=0)

    # convert the features and targets to the 32-bit format suitable for the
    # model
    x_data = x_data.astype(np.float32)
    y_data = y_data.astype(np.int32)

    # Split data into training and classifying:
    ind = int(0.8 * len(x_data))
    exp = train(x_data[:ind], y_data[:ind])

    for output in classify(exp, x_data[ind:], y_data[ind:]):
        print output


if __name__ == '__main__':
    main()
