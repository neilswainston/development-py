'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import json

from synbiochem.utils import seq_utils


query = 'hexokinase and taxonomy:9606'

fields = ['entry name', 'protein names', 'sequence', 'ec', 'organism',
          'organism-id', 'database(PDB)']

result = seq_utils.search_uniprot(query, fields)

print json.dumps(result, indent=4, sort_keys=True)
