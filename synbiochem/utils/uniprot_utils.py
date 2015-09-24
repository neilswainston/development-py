'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import csv
import urllib2


def get_sequences(uniprot_ids):
    '''Gets dictionary of ids to sequences from Uniprot.'''
    query = '+or+'.join(['id:' + uniprot_id for uniprot_id in uniprot_ids])
    url = 'http://www.uniprot.org/uniprot/?query=' + query + \
        '&format=tab&columns=id,sequence'
    return dict(list(csv.reader(urllib2.urlopen(url), delimiter='\t'))[1:])
