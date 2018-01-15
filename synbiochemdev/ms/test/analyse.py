'''
synbiochem (c) University of Manchester 2017

All rights reserved.

@author: neilswainston
'''
# pylint: disable=invalid-name
import re
import sys

import pandas as pd


def analyse(df):
    '''analyse.'''
    result_df = df.groupby(['plasmid',
                            'strain',
                            'treatment'])['mg/L'].agg({'mg/L': ['mean',
                                                                'std']})

    print result_df


def get_data(filename):
    '''get data.'''
    df = pd.read_csv(filename)

    df['Name'] = df['Name'].apply(_fix_name)
    names_df = df['Name'].str.split('[ _]', expand=True)
    names_df.columns = ['strain', 'plasmid', 'mutant', 'colony', 'media',
                        'carbon source', 'antibiotic', 'treatment',
                        'injection']
    names_df['strain'] = names_df['strain'] + '_' + names_df['mutant']
    names_df.index = df.index

    df = df.join(names_df)

    df.to_csv('out.csv')

    return df


def _fix_name(name):
    '''fix name.'''
    if name.count('_') == 0 or name.count('_') == 6:
        return name

    pos = [m.start() for m in re.finditer('_', name)][1]
    return name[:pos] + '_' + name[pos:]


def main(args):
    '''main method.'''
    df = get_data(args[0])
    analyse(df)


if __name__ == '__main__':
    main(sys.argv[1:])
