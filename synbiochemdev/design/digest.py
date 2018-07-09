'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys

from Bio.Restriction.Restriction import RestrictionBatch
from Bio.Seq import Seq
from synbiochem.utils import ice_utils


def digest(ice_url, ice_username, ice_password, ice_ids):
    '''Digest sequences from ICE.'''
    ice_client = ice_utils.ICEClient(ice_url, ice_username, ice_password)
    rest_batch = RestrictionBatch(['BamHI', 'EcoRI'])

    for ice_id in ice_ids:
        seq = ice_client.get_ice_entry(ice_id).get_seq()
        rest_sites = rest_batch.search(Seq(seq), linear=False)
        rest_sites = sorted([0] + [site - 1
                                   for sites in rest_sites.values()
                                   for site in sites])

        dig_seqs = [seq[i:j]
                    for i, j in zip(rest_sites, rest_sites[1:] + [None])]

        dig_seqs = dig_seqs[1:-1] + [dig_seqs[-1] + dig_seqs[0]]

        print '\t'.join([ice_id] + [str(len(dig_seq))
                                    for dig_seq in dig_seqs])


def main(args):
    '''main method.'''
    digest(args[0], args[1], args[2], args[3:])


if __name__ == '__main__':
    main(sys.argv[1:])
