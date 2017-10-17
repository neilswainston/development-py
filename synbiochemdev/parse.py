'''
(c) GeneGenie 2017

All rights reserved.

@author: neilswainston
'''
import sys

import pandas as pd


def main(args):
    '''main method.'''
    df_format = pd.read_csv(args[0], sep='\t')
    df_format = df_format[df_format.Description != 'Unused']

    df_data = pd.DataFrame(columns=df_format['Field Name'])

    with open(args[1]) as datafile:
        for line in datafile:
            row = {field: _get_val(line, col, lngth, typ)
                   for field, col, lngth, typ in zip(df_format['Field Name'],
                                                     df_format['Column'],
                                                     df_format['Length'],
                                                     df_format['Type'])}
            df_data = df_data.append(pd.Series(row), ignore_index=True)

    df_data.to_csv(args[2], sep='\t')


def _get_val(line, col, lngth, typ):
    '''Get value.'''
    return line[col - 1:col - 1 + lngth]


if __name__ == '__main__':
    main(sys.argv[1:])
