'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name

import io
import sys

from synbiochem.utils import ice_utils, dna_utils

import pandas as pd


def write_plasmid(filename, ice_url, ice_username, ice_password,
                  group_name=None):
    '''Writes plasmids.'''
    df = pd.read_csv(filename)
    ice_client = ice_utils.ICEClient(ice_url, ice_username, ice_password)
    plasmid_ids = []

    for _, row in df.iterrows():
        part = ice_client.get_ice_entry(row['part'])
        vector = ice_client.get_ice_entry(row['vector'])

        # Form name:
        name = part.get_metadata()['name'] + \
            ' (' + vector.get_metadata()['name'] + ')'

        # Create plasmid:
        plasmid = ice_utils.ICEEntry(typ='PLASMID')
        plasmid.set_values({'name': name[:127], 'shortDescription': name})
        plasmid.set_dna(dna_utils.concat([part.get_dna(), vector.get_dna()]))

        # Upload to ICE:
        ice_client.set_ice_entry(plasmid)

        # Add links:
        ice_client.add_link(plasmid.get_ice_id(), part.get_ice_id())
        ice_client.add_link(plasmid.get_ice_id(), vector.get_ice_id())
        ice_client.set_ice_entry(plasmid)

        # Add groups:
        if group_name:
            groups = ice_client.get_groups()
            ice_client.add_permission(plasmid.get_ice_id(),
                                      groups[group_name])

        plasmid_ids.append(plasmid.get_ice_id())

    # Update dataframe:
    df['plasmid'] = plasmid_ids

    df.to_csv('out.csv', index=False)


def set_groups(filename, ice_url, ice_username, ice_password, group_name):
    '''Sets group permissions.'''
    with io.open(filename, 'r', newline='\r') as fle:
        ice_client = ice_utils.ICEClient(ice_url, ice_username, ice_password)
        groups = ice_client.get_groups()

        for ice_id in fle:
            ice_client.add_permission(ice_id, groups[group_name])


def main(args):
    '''main method.'''
    write_plasmid(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
