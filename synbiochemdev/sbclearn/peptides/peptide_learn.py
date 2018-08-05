'''
sbclearn (c) University of Manchester 2017

sbclearn is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=no-member
# pylint: disable=ungrouped-imports
import sys

from sklearn import preprocessing

import numpy as np
import pandas as pd
from sbclearn.utils.validator import k_fold_cross_valid


def get_data(filename):
    '''Gets data.'''
    df = pd.read_table(filename)

    # Liquid required at pH 4:
    df.loc[df['Physical state pH4'] == 'L', 'Physical state pH4'] = 0.0
    df.loc[df['Physical state pH4'] == 'V', 'Physical state pH4'] = 0.5
    df.loc[df['Physical state pH4'] == 'G', 'Physical state pH4'] = 1.0

    # Gel required at pH 7:
    df.loc[df['Physical state pH7'] == 'G', 'Physical state pH7'] = 0.0
    df.loc[df['Physical state pH7'] == 'V', 'Physical state pH7'] = 0.5
    df.loc[df['Physical state pH7'] == 'L', 'Physical state pH7'] = 0.75
    df.loc[df['Physical state pH7'] == 'P', 'Physical state pH7'] = 1.0

    # Normalise:
    num = preprocessing.MinMaxScaler(
        feature_range=(0.1, 0.9)).fit_transform(df.ix[:, 2:])

    # Maximise:
    df_ser = pd.DataFrame(num).applymap(lambda x: 1 - x)

    # Reform DataFrame:
    df_ser.columns = df.ix[:, 2:].columns
    df = pd.concat([df.ix[:, :2], df_ser], axis=1)

    # Set objective:
    df['obj'] = df.ix[:, 2:].prod(axis=1)

    x_data = np.array(sbclearn.get_aa_props(df['Sequence'].tolist()))
    y_data = df['obj']

    return x_data, y_data


def main(args):
    '''main method.'''
    x_data, y_data = get_data(args[0])
    results = k_fold_cross_valid((x_data, y_data))
    sbclearn.plot(results, 'Prediction of peptide fitness')


if __name__ == '__main__':
    main(sys.argv[1:])
