'''
Created on 29 Nov 2016

@author: neilswainston
'''
# pylint: disable=invalid-name
import sys

from synbiochem.utils import ice_utils

import pandas as pd


def write(filename, ice_url, ice_username, ice_password,
          group_name=None):
    '''Writes strains.'''
    ice_client = ice_utils.ICEClient(ice_url, ice_username, ice_password)
    df = pd.read_csv(filename)

    for row in df.T.to_dict().values():
        entry = ice_utils.ICEEntry(typ='PART')
        entry.set_parameter('Type', row.pop('Type'))
        entry.set_values(row)
        ice_id = ice_client.set_ice_entry(entry)

        if group_name:
            groups = ice_client.get_groups()
            ice_client.add_permission(ice_id, groups[group_name])

        print ice_id


def main(args):
    '''main method.'''
    write(*args)


if __name__ == '__main__':
    main(sys.argv[1:])
