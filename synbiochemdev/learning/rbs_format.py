'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import collections
from synbiochem.utils import sequence_utils

with open('rbs.txt', 'rU') as infile:
    id_seqs = {idx: line.split('\t')[0]
               for idx, line in enumerate(infile)}
    id_seqs = collections.OrderedDict(sorted(id_seqs.items()))

    sequence_utils.write_fasta(id_seqs, 'rbs_out.txt')
