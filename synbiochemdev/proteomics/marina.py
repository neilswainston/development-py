'''
Created on 5 Jan 2018

@author: neilswainston
'''
# pylint: disable=invalid-name
import sys

from synbiochem.utils import seq_utils

import pandas as pd


def main(args):
    '''main method.'''
    df = pd.read_csv(args[0])

    df['Symbol'] = df['Symbol'].str.split().str[0]
    df['Symbol'].fillna('', inplace=True)
    df['uniprot'] = df['Symbol'].map(get_uniprot_id)

    df.to_csv('out.csv', index=False)


def get_uniprot_id(gene_id):
    '''Get uniprot id.'''
    result = seq_utils.search_uniprot(
        'gene:' + gene_id + '+AND+organism:9606', [], 1)

    if result:
        print result
        return result[0]['Entry']

    return None


if __name__ == '__main__':
    main(sys.argv[1:])
