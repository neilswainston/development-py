#  -*- coding: utf-8 -*-
'''
Spyder Editor

This is a temporary script file.
'''

# import climate # for automatic logging
from compiler.ast import flatten
import os
import time

from sklearn.metrics import log_loss
from sklearn.metrics import roc_auc_score as auc
import theanets

import pandas


# from sklearn.metrics import confusion_matrix
# import matplotlib.pyplot as plt
# This using input as numerical codes, 0,1,2,3... should really use
# 'One-hot-vectors' - this would be done using pandas.get_dummies
# See:
# http://pandas.pydata.org/pandas-docs/version/0.17.0/generated/pandas.get_dummies.html
def get_data(f_name, cat_cols, excluded_cols, target):
    '''Read data from csv, and use the 'target' columns as y-data,
    everything else is x-data.'''

    # perhaps we need to return the codes too so we can go from
    # numbers back to categories...
    raw_data = pandas.read_csv(f_name)
    x_data = raw_data.drop(excluded_cols, 1)

    # Convert strings to categories
    for col in cat_cols:
        x_data[col] = x_data[col].astype('category')
        x_data[col] = x_data[col].cat.codes

    y_data = raw_data[target]
    y_data = y_data.astype('category')
    y_data = y_data.cat.codes
    y_data = y_data.astype('int32')

    # Theanets/GPU want's imput x's as float32, y's as int32
    x_data = x_data.astype('float32')
    return (x_data, y_data)


def do_net_run(train, valid, test, act_type, num_hidden, num_nodes, num_in,
               num_out, drop_out=0.0, num_reps=1):
    '''Perform classification.'''
    layer_def = dict(size=num_nodes, activation=act_type)

    net_def = flatten(
        [num_in, [layer_def for _ in range(num_hidden)], num_out])

    net = theanets.Classifier(layers=net_def)
    net.save('model1')

    for rep in range(num_reps):
        net = theanets.Classifier.load('model1')
        start = time.clock()
        epoch = 0
        for trn, vdl in net.itertrain(train, valid,
                                      # stopping criteria
                                      # how many 'fails' we allow.
                                      patience=5,
                                      # validatations must improve by x%
                                      min_improvement=0.005,
                                      # Do validation every 'n' steps
                                      validate_every=1,
                                      # batch learning
                                      batch_size=5,
                                      # minibatches per epoch
                                      # train_batches=30,
                                      # dropouts etc
                                      hidden_dropout=drop_out,
                                      input_dropout=drop_out,
                                      algo='rmsprop'):
            # print 'Tr loss: ' + str(trn['loss']), 'Vld loss: ',
            # str(vdl['loss'])
            epoch = epoch + 1

        duration = time.clock() - start

        x_data, y_data = test
        auroc = auc(y_data, net.predict(x_data))
        tst_loss = log_loss(y_data, net.predict_proba(x_data))

        x_data, y_data = train
        trn_loss = log_loss(y_data, net.predict_proba(x_data))

        x_data, y_data = valid
        vld_loss = log_loss(y_data, net.predict_proba(x_data))

        print '\t'.join([str(x) for x in [act_type, num_hidden, num_nodes,
                                          drop_out, rep, epoch, trn_loss,
                                          vld_loss, tst_loss, auroc,
                                          duration]])


def get_all_data():
    '''Gets all data.'''
    #  use unix directory '/' seperator to avoid python
    #  taking '\' to be the escape character
    # data_path='d:/data/UCI_ML_data/'
    data_path = os.path.dirname(os.path.realpath(__file__))

    # list of catagorical column -> to be converted to numeric codes.
    catcol = ['cap-shape', 'cap-surface', 'cap-color', 'bruises', 'odor',
              'gill-attachment', 'gill-spacing', 'gill-size', 'gill-color',
              'stalk-shape', 'stalk-root', 'stalk-surface-above-ring',
              'stalk surface-below ring', 'stalk-color-above-ring',
              'stalk-color-below-ring', 'veil-type', 'veil-color',
              'ring-number', 'ring-type', 'spore-print-color', 'population',
              'habitat']

    # list of columns to excluded (not 'x' features)
    ex_cols = ['E/P']
    # target column i.e. 'y' - assumed categorical
    target = 'E/P'

    train = get_data(
        os.path.join(data_path, 'mush_train20.csv'), catcol, ex_cols, target)
    valid = get_data(
        os.path.join(data_path, 'mush_valid20.csv'), catcol, ex_cols, target)
    test = get_data(
        os.path.join(data_path, 'mush_test20.csv'), catcol, ex_cols, target)

    return train, valid, test


def scan(act_list, hid_lay_list, neuron_count_list, dropout_list, reps):
    '''Perform parameter scan.'''
    train, valid, test = get_all_data()

    print '\t'.join(['act_type', 'num_hidden', 'num_nodes', 'drop_out', 'Rep',
                     'Numepochs', 'TrLoss', 'VldLoss', 'TestLoss', 'AUR',
                     'Time'])

    for actv in act_list:
        for hid_lay in hid_lay_list:
            for neuron_count in neuron_count_list:
                for drp in dropout_list:
                    do_net_run(train, valid, test, actv, hid_lay, neuron_count,
                               22, 2, drp, reps)


def main():
    '''Main method.'''
    act_list = ['relu', 'prelu', 'lgrelu']
    hid_lay_list = [1, 2, 3, 5, 10]
    neuron_count_list = [20, 50]
    dropout_list = [0, 0.1, 0.3, 0.5, 0.7]
    reps = 3

    #  something on pre-normalisation

    scan(act_list, hid_lay_list, neuron_count_list, dropout_list, reps)


if __name__ == '__main__':
    main()
