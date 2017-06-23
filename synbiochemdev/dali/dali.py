'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import csv
import re
import sys
import tempfile
import urllib

from Bio.PDB import parse_pdb_header
import requests

import pandas as pd


def analyse(filename):
    '''Analyse method.'''
    df = pd.read_table(filename)
    cache = {}
    pdb_data = []

    for chain in df['Chain']:
        pdb_id = chain.split('-')[0]
        fields = cache.get(pdb_id, None)

        if not fields:
            url = 'https://files.rcsb.org/download/%s.pdb' % pdb_id
            tmpfile = tempfile.NamedTemporaryFile(delete='False')
            urllib.urlretrieve(url, tmpfile.name)

            with open(tmpfile.name, 'r') as fle:
                header_dict = parse_pdb_header(fle)
                fields = header_dict['source']['1']
                fields['uniprot'] = _get_uniprot_id(pdb_id)

                if fields['uniprot']:
                    fields['uniprot url'] = '=HYPERLINK("' + \
                        'http://www.uniprot.org/uniprot/' + \
                        fields['uniprot'] + '")'

                fields['ec_number'] = \
                    header_dict['compound']['1'].get('ec_number', '')
                fields['journal'] = header_dict['journal']

                dois = re.findall(r'10.\d{4,9}/[-._;()/:A-Z0-9]+',
                                  fields['journal'])

                if dois:
                    fields['doi url'] = '=HYPERLINK("' + \
                        'https://dx.doi.org/' + dois[0] + '")'

                cache[pdb_id] = fields

        pdb_data.append(fields)

    pdb_df = pd.DataFrame(pdb_data)

    df = pd.concat([df, pdb_df], axis=1)
    df.to_csv(filename + '_out.xls', sep='\t', quoting=csv.QUOTE_NONE)


def _get_uniprot_id(pdb_id):
    '''Gets uniprot ids.'''
    url = 'http://www.uniprot.org/uniprot/?query=database:(type:pdb ' + \
        pdb_id + ')&sort=score&columns=id&format=tab'

    resp = requests.get(url)

    if resp.status_code is not 200:
        resp.raise_for_status()

    for row in csv.reader(resp.text.splitlines()[1:],
                          delimiter='\t'):
        return row[0]

    return None


def main(args):
    '''main method.'''
    for filename in args:
        analyse(filename)


if __name__ == '__main__':
    main(sys.argv[1:])
