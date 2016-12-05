'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys

from synbiochem.utils.ice_utils import ICEClient
from synbiochem.utils.net_utils import NetworkError


def analyse(url, username, password):
    '''Exports a plasmids constituent parts for ordering.'''
    ice_client = ICEClient(url, username, password)
    deleted = 0

    for ice_number in range(0, 2358):
        try:
            ice_entry = ice_client.get_ice_entry(ice_number)
            metadata = ice_entry.get_metadata()

            if metadata['visible'] == 'OK':
                seq = ice_entry.get_seq()
                results = [metadata['partId'],
                           metadata['name'] if 'name' in metadata else '',
                           metadata['type'] if 'type' in metadata else '',
                           metadata['visible']
                           if 'visible' in metadata else '',
                           ice_entry.get_parameter('Type'),
                           len(seq) if seq is not None else 0]

                print '\t'.join([str(res) for res in results])
        except NetworkError:
            # Take no action
            deleted += 1
        except UnicodeEncodeError, err:
            print err


def main(args):
    '''main method.'''
    analyse(args[0], args[1], args[2])


if __name__ == '__main__':
    main(sys.argv[1:])
