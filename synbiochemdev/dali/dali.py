'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import csv
import sys
import tempfile
import urllib

from Bio.PDB import parse_pdb_header

import pandas as pd


def analyse(filename):
    '''Analyse method.'''
    df = pd.read_table(filename)
    pdb_data = []

    count = 0

    for chain in df['Chain']:
        count = count + 1

        if count == 5:
            break

        pdb_id = chain.split('-')[0]
        url = 'https://files.rcsb.org/download/%s.pdb' % pdb_id
        tmpfile = tempfile.NamedTemporaryFile(delete='False')
        urllib.urlretrieve(url, tmpfile.name)

        with open(tmpfile.name, 'r') as fle:
            header_dict = parse_pdb_header(fle)
            fields = header_dict['source']['1']
            fields['ec_number'] = header_dict['compound']['1'].get('ec_number',
                                                                   '')
            fields['journal'] = header_dict['journal']
            pdb_data.append(fields)

    df = pd.concat([df, pd.DataFrame(pdb_data)], axis=1)
    df.to_csv(filename + '_out.txt', sep='\t', quoting=csv.QUOTE_NONE)


def main(args):
    '''main method.'''
    for filename in args:
        analyse(filename)


if __name__ == '__main__':
    main(sys.argv[1:])
