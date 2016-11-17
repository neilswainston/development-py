# -*- coding: utf-8 -*-
'''
Created on Wed Jan 13 11:54:49 2016

@author: mcdssso
'''


# Using holygrain Jan 2016 version, also needs synbiochem
# https://github.com/neilswainston/HolyGrail/tree/master/holygrail
# https://github.com/synbiochem/synbiochem-py
import inspect
import itertools
import os.path
import random
import time

from holygrail import get_input_data
from sklearn.metrics import accuracy_score as accs
from sklearn.metrics import confusion_matrix as confmat
from sklearn.metrics import log_loss
import cPickle as pkl
import holygrail.data as hgd
import numpy as np
import pandas as pd
import theanets as thts


def getData(numSamples, rgx, rgLbl, min_ham=3, trn=2, val=1, tst=1):
    '''Retreives peptide string matching regular expression list and randomly
    samples them
    into training, validation and test sets.
    :param numSamples: number of samples per category.
    :param rgx: regular expression matching peptide secondary structures to
    retrieve.
    :param rgLbl: regex label, replacement labels used for column headers
    :param min_ham: minimum hamming distance between any two peptides.
    :param trn,val,tst: correspond to training:validation:test ratio
    :returns: tupple, (train,valid,test), each element is a pandas DataFrame.
    '''

    sq, _ = hgd.sample_seqs(numSamples, rgx, min_hamming=min_ham)

    dt = pd.DataFrame()
    for key in sq.iterkeys():
        tmp = pd.DataFrame(sq[key])
        tmp['Class'] = key
        dt = pd.concat([dt, tmp], ignore_index=True)

    dt.columns = ['seq', 'struct', 'PepChainID', 'StartStop', 'Class']

    dt[['PepID', 'ChainID']] = dt['PepChainID'].apply(pd.Series)
    dt[['pStart', 'pStop']] = dt['StartStop'].apply(pd.Series)

    dt.drop(['PepChainID', 'StartStop'], axis=1, inplace=True)

    # humanise class names in Class column
    for idx, key in enumerate(rgx):
        dt.loc[dt.Class == key, 'Class'] = rgLbl[idx]

    dt = dt.iloc[np.random.permutation(len(dt))]

    dt.reset_index(drop=True, inplace=True)

    tot = trn + val + tst

    rows, _ = np.shape(dt)

    k1 = int(rows * trn / tot)
    k2 = int(rows * (trn + val) / tot)

    trDat = dt.iloc[0:k1]
    vlDat = dt.iloc[k1:k2]
    tstDat = dt.iloc[k2:]

    return trDat, vlDat, tstDat


def getRndTstRow(dat):
    x, y = dat
    nRows, inLen = np.shape(x)
    idx = np.random.randint(0, nRows)
    return np.reshape(x[idx, :], (1, inLen)), y[idx]


def mungeData4NN(dat):
    '''Convert peptide data into format suitable neural network.
    Peptides are represented as the five amino-acid properties KD
    Hydrophobicity,
    EIIP (electron-ion interaction ptential), Helix, Sheet, and Turn for each
    residue.

    :param dat: Pandas Dataframe containing peptide data
    :returns: (xdat,ydata) numpy arrays suitable for input into a Theanets NN
    running on GPU.'''

    # cut out un-needed columns
    tmp = dat.loc[:, ['seq', 'struct', 'Class']]

    rowCount, _ = np.shape(tmp)

    # TODO: This could be generalised to cope with peptides of different
    # lengths.
    # If peptides are not same length, one would pad to the max length with
    # blancks,
    # One could also add an end marker stop/start character.
    pepLength = 9
    inputLength = pepLength * 5

    xs = np.zeros((rowCount, inputLength), 'f')

    j = 0
    for _, row in tmp.iterrows():
        xs[j, :] = np.reshape(get_input_data(row.seq), (1, 45))
        j += 1
    ys = pd.Categorical(tmp.Class, categories=regxLabels).codes.astype('int32')
    return xs, ys


def doNetRun(DoPretrain, actType, numHidden, numNodes, dropOut, NumReps,
             AdjustWforDropout, L1, L2,
             LearningRate, Momentum, Algorithm, maxUpdate, BatchSize, Patience,
             MinImprovement, ValidateEveryN):
    '''Performs *numRep* neural network runs, saving the best run so far
    (global var)
    :params: NN hyper-parameters.
    :returns: run parameters + performance measures (DataFrame).'''

    layerDef = dict(size=numNodes, activation=actType)

    netDef = []
    netDef.append(numFeatures)
    for _ in range(numHidden):
        netDef.append(layerDef)
    netDef.append(numCats)

    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    saveCols = []
    for arg in args:
        print arg, ': ', values[arg], ',',
        saveCols.append(arg)
    print

    saveCols = saveCols + \
        ['TrLoss', 'VldLoss', 'TstLoss', 'TestAcc', 'Time', 'BestEpoch']

    hist = pd.DataFrame(columns=saveCols)

    for _ in range(NumReps):
        global kountRuns
        kountRuns = kountRuns + 1

        # use new seed for each run
        ii32 = np.iinfo(np.int32)
        seed = random.randint(0, ii32.max)
        print 'seed: ', seed
        # although numpy.random.seed should be uint32, Theanets only checks for
        # int
        # and fails if given a uint32 stored as a long in Python

        net = thts.Classifier(layers=netDef, rng=seed)

        t0 = time.clock()
        Epoch = 0
        if DoPretrain:
            print('Train phase I:')
            net.train(train, valid,
                      patience=Patience,
                      learning_rate=LearningRate,
                      momentum=Momentum,
                      min_improvement=MinImprovement,
                      validate_every=ValidateEveryN,
                      max_updates=maxUpdate,
                      input_dropout=dropOut,
                      hidden_dropout=dropOut,
                      algo='layerwise',
                      weight_l1=L1,  # L1 norm sparsity
                      weight_l2=L2,  # L2 norm weight decay
                      batch_size=BatchSize)

        print('Train phase II:')
        Epoch = 0
        lastLoss = np.Inf
        lastEpoch = 0
        for tr, vl in net.itertrain(train, valid,
                                    patience=Patience,
                                    learning_rate=LearningRate,
                                    momentum=Momentum,
                                    min_improvement=MinImprovement,
                                    validate_every=ValidateEveryN,
                                    max_updates=maxUpdate,
                                    input_dropout=dropOut,
                                    hidden_dropout=dropOut,
                                    algo=Algorithm,
                                    weight_l1=L1,  # L1 norm sparsity
                                    weight_l2=L2,  # L2 norm weight decay
                                    batch_size=BatchSize):
            Epoch = Epoch + 1
            vloss = vl['loss']

            if (lastLoss - vloss) >= MinImprovement:
                lastLoss = vloss
                lastEpoch = Epoch
                flg = ' *' + str(lastEpoch)
            else:
                flg = ''
            print Epoch, 'trLoss: %.4f' % tr['loss'], ' vlLoss: %.4f' % vloss, \
                ' vlacc: %.4f' % vl['acc'], flg

        t1 = time.clock() - t0
        print 'Time: ', t1, ' Epochs:', Epoch

        if AdjustWforDropout:
            fact = 1.0 - dropOut
            for ll in net.layers:
                if (ll.name != 'in'):
                    w = net.find(ll.name, 'w')
                    w.set_value(w.get_value() * fact)

        X, y = train
        trnLoss = log_loss(y, net.predict_proba(X))

        X, y = valid
        vldLoss = log_loss(y, net.predict_proba(X))

        X, y = test
        ypp = net.predict_proba(X)
        yp = net.predict(X)
        acc = net.score(X, y)
        tstLoss = log_loss(y, ypp)

        print Epoch, 'trLoss: %.4f' % trnLoss, ' vlLoss: %.4f' % vldLoss, \
            ' vlacc: %.4f' % acc
        print 'Best Epoch: ', lastEpoch

        cf = confmat(y, yp)
        print 'Test-set confusion matrix:'
        print cf

        global bestLoss
        global bestParams

        dta = dict()

        for arg in args:
            dta[arg] = values[arg]

        dta['TrLoss'] = trnLoss
        dta['VldLoss'] = vldLoss
        dta['TstLoss'] = tstLoss
        dta['TestAcc'] = acc
        dta['Time'] = t1
        dta['BestEpoch'] = lastEpoch

        nr = pd.DataFrame([dta])

        if (tstLoss <= bestLoss):
            bestLoss = tstLoss
            net.save('bestModel')
            bestParams = nr

        hist = hist.append(nr, ignore_index=True)

        # re-order columns...
        hist = hist[saveCols]

    return hist


def doBigLoops(JustTest=True):
    '''Do the NN main loop.'''
    global bestLoss
    global bestModel

    bestLoss = np.Inf
    bestModel = None

    if JustTest:
        pretrainList = [False]
        # http://theanets.readthedocs.org/en/stable/api/activations.html
        actList = ['relu']
        HLList = [2]
        neuronCountList = [400]
        dropoutList = [0.0]
        repsLst = [1]
        # this actually does weight adjustment before predict
        doDropOutOnPredictLst = [False]
        L1WList = [0.0]
        L2WList = [0.0]
        lRateList = [0.005]
        momentumList = [0.9]
        algoList = ['nag']  # sgd,nag,rprop,rmsprop,adadelta,esgd,adam
        maxUpdLst = [300]
        BatchSizeLst = [64]
        PatienceLst = [5]
        min_impLst = [0.002]
        ValidateEveryNLst = [1]
    else:
        pretrainList = [False]
        actList = ['relu']
        HLList = [2]
        neuronCountList = [1000]
        dropoutList = [0.0]
        repsLst = [1]
        # this actually does weight adjustment before predict
        doDropOutOnPredictLst = [False]
        L1WList = [0.0]
        L2WList = [0.0]
        lRateList = [0.001]
        momentumList = [0.7]
        algoList = ['nag']
        maxUpdLst = [300]
        BatchSizeLst = [256]
        PatienceLst = [5]
        min_impLst = [0.005]
        ValidateEveryNLst = [1]

    zTime = time.clock()

    runsToDo = 0

    for args in itertools.product(pretrainList, actList, HLList,
                                  neuronCountList,
                                  dropoutList, repsLst, doDropOutOnPredictLst,
                                  L1WList, L2WList,
                                  lRateList, momentumList, algoList, maxUpdLst,
                                  BatchSizeLst,
                                  PatienceLst, min_impLst, ValidateEveryNLst):
        runsToDo += 1

    global kountRuns

    kountRuns = 1

    resultData = None

    for args in itertools.product(pretrainList, actList, HLList,
                                  neuronCountList,
                                  dropoutList, repsLst, doDropOutOnPredictLst,
                                  L1WList, L2WList,
                                  lRateList, momentumList, algoList, maxUpdLst,
                                  BatchSizeLst,
                                  PatienceLst, min_impLst, ValidateEveryNLst):
        # print args
        print 'Run: ', str(kountRuns), ' of ', str(runsToDo)
        tmp = doNetRun(*args)
        if resultData is not None:
            resultData = resultData.append(tmp)
        else:
            resultData = tmp
        print 'Time elapsed ', (time.clock() - zTime)

    return resultData


def makeBestPredictionTable(dat, label):
    net = thts.Classifier.load('bestModel')
    X, y = dat
    yp = net.predict(X)
    acc = accs(y, yp)
    tstLoss = log_loss(y, net.predict_proba(X))

    print label, ' Loss: %.4f' % tstLoss, ' Acc: %.4f' % acc

    dtst = pd.DataFrame(X)
    dtst.columns = ['x' + str(i + 1) for i in range(45)]
    dtst['Type'] = label
    dtst['Acc'] = acc
    dtst['loss'] = tstLoss
    dtst['y'] = y
    dtst['yp'] = yp
    return dtst


aaCodes = list('ABCDEFGHIKLMNPQRSTUVWXYZ ')
ssCodes = list('HBEGITS ')

# regex and label should in same order!
regxStrings = ['.{4}[HIG].{4}', '.{4}[EB].{4}', '.{4}[TSL ].{4}']
regxLabels = ['HIG', 'EB', 'TSL']

# regxStrings=['H{9}','E{9}','[^EH]{9}']
# regxLabels=['H','E','O']

samplesPerCat = 20000

fileName = 'pep60K.p'
# fileName='pepHEO60k.p'

forceRegen = False

doSave = True
if (os.path.exists(fileName) and not forceRegen):
    print 'Loading data.'
    doSave = False
    tr, vl, ts = pkl.load(open(fileName, 'rb'))
else:
    print 'Generating data... may take some time.'
    tmp = getData(samplesPerCat, regxStrings, regxLabels, min_ham=3)
    pkl.dump(tmp, open(fileName, 'wb'))
    tr, vl, ts = tmp
    del tmp

# Save combined data (if new)
if doSave:
    t1 = tr
    t2 = vl
    t3 = ts

    t1['Type'] = 'train'
    t2['Type'] = 'valid'
    t3['Type'] = 'test'

    dt = pd.concat([t1, t2, t3], ignore_index=True)
    del t1, t2, t3

    dt['nClass'] = pd.Categorical(
        dt['Class'], categories=regxLabels).codes.astype('int32')

    dt[dt.Type == 'test'].to_csv('just_test.csv')

    dt.to_csv(fileName + '.combinedData.csv')
    del dt

train = mungeData4NN(tr)
valid = mungeData4NN(vl)
test = mungeData4NN(ts)
del ts, vl, tr

rowcount, numFeatures = train[0].shape

numCats = len(np.unique(train[1]))


print 'Starting loops...\n'
resultData = doBigLoops(JustTest=True)
resultData.to_csv(fileName + '.ResultSet.csv', index=False)


# Save 'best' predictions.
t1 = makeBestPredictionTable(train, 'train')
t2 = makeBestPredictionTable(valid, 'valid')
t3 = makeBestPredictionTable(test, 'test')

RDF = pd.concat([t1, t2, t3])
del t1, t2, t3

RDF.to_csv(fileName + '.BestResult.csv')

print 'Run complete.'
