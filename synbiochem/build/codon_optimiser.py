'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import operator
import random
import re
import urllib2
import synbiochem.utils


class CodonOptimiser(object):
    '''Class to support codon optimisation.'''

    def __init__(self, taxonomy_id):
        self.__taxonomy_id = taxonomy_id
        self.__codon_usage_table = self.__get_codon_usage_table()

    def optimise(self, protein_seqs):
        '''Codon optimises the supplied protein sequences.'''
        optimised_seqs = []

        for protein_seq in protein_seqs:
            optimised_seqs.append(''.join([self.get_random_codon(aa)
                                           for aa in protein_seq]))

        return optimised_seqs

    def get_random_codon(self, amino_acid):
        '''Returns a random codon for a given amino acid,
        based on codon probability from the codon usage table.'''
        codon_usage = self.__codon_usage_table[amino_acid]
        rand = random.random()
        cumulative_prob = 0

        for codon, prob in iter(reversed(codon_usage)):
            cumulative_prob += prob

            if cumulative_prob > rand:
                return codon

    def __get_codon_usage_table(self):
        '''Gets the codon usage table for a given taxonomy id.'''
        codon_usage_table = {aa_code: {}
                             for aa_code in synbiochem.utils.AA_CODES.values()}

        url = 'http://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=' \
            + self.__taxonomy_id + '&aa=1&style=GCG'

        in_codons = False

        for line in urllib2.urlopen(url):
            if line == '<PRE>\n':
                in_codons = True
            elif line == '</PRE>\n':
                break
            elif in_codons:
                values = re.split('\\s+', line)

                if values[0] in synbiochem.utils.AA_CODES:
                    aa_code = synbiochem.utils.AA_CODES[values[0]]
                    codon_usage = codon_usage_table[aa_code]
                    codon_usage[values[1]] = float(values[3])

        codon_usage_table.update((x, _scale(y))
                                 for x, y in codon_usage_table.items())

        print codon_usage_table

        return codon_usage_table


def _scale(codon_usage):
    '''Scale codon usage values to add to 1.'''
    codon_usage = dict([(key, value / sum(codon_usage.values()))
                        for key, value in codon_usage.items()])

    return sorted(codon_usage.items(), key=operator.itemgetter(1))
