'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author: neilswainston
'''
import sys

from synbiochem.utils import ice_utils, sbol_utils


class SequenceAnalyser(object):
    '''Class to analyse sequences.'''

    def __init__(self, url, username, psswrd):
        self.__ice_client = ice_utils.ICEClient(url, username, psswrd)

    def digest(self, ice_id, restricts, circular=False):
        '''Performs digest on supplied ICE entry id.'''
        ice_entry = self.__ice_client.get_ice_entry(ice_id)
        return [sbol_utils.get_seq(doc)
                for doc in sbol_utils.apply_restricts(ice_entry.get_sbol_doc(),
                                                      restricts,
                                                      circular=circular)]

    def pcr(self, ice_id, for_primer, rev_primer, circular=False):
        '''Performs digest on supplied ICE entry id.'''
        ice_entry = self.__ice_client.get_ice_entry(ice_id)
        return [sbol_utils.get_seq(doc)
                for doc in sbol_utils.apply_pcr(ice_entry.get_sbol_doc(),
                                                for_primer, rev_primer,
                                                circular=circular)]


def main(argv):
    '''Main method.'''
    seq_anal = SequenceAnalyser(argv[0], argv[1], argv[2])

    for ice_id in argv[3:-2]:
        digest = [len(seq) for seq in seq_anal.digest(ice_id,
                                                      argv[-2].split(','),
                                                      circular=True)]
        pcr = [len(seq) for seq in seq_anal.pcr(ice_id,
                                                *argv[-1].split(','),
                                                circular=True)]
        print '\t'.join([ice_id, str(digest), str(pcr)])


if __name__ == '__main__':
    main(sys.argv[1:])
