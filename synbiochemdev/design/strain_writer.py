'''
Created on 29 Nov 2016

@author: neilswainston
'''
import io
import sys

from synbiochem.utils import ice_utils


def write_strain(filename, ice_url, ice_username, ice_password,
                 group_name=None):
    '''Writes strains.'''
    with io.open(filename, 'r', newline='\r') as fle:
        pairs = [str(line).split() for line in fle]

        ice_client = ice_utils.ICEClient(ice_url, ice_username, ice_password)

        for pair in pairs:
            host = ice_client.get_ice_entry(pair[0])
            plasmid = ice_client.get_ice_entry(pair[1])
            host_metadata = host.get_metadata()
            plasmid_metadata = plasmid.get_metadata()
            name = host_metadata['name'] + \
                ' (' + plasmid_metadata['name'] + ')'

            strain = ice_utils.ICEEntry(typ='STRAIN')
            strain.set_values({'name': name, 'shortDescription': name})
            strain.set_parameter('Taxonomy', host.get_parameter('Taxonomy'))
            ice_client.set_ice_entry(strain)

            ice_client.add_link(strain.get_ice_id(), host.get_ice_id())
            ice_client.add_link(strain.get_ice_id(), plasmid.get_ice_id())
            ice_client.set_ice_entry(strain)

            if group_name:
                groups = ice_client.get_groups()
                ice_client.add_permission(strain.get_ice_id(),
                                          groups[group_name])

            print host.get_ice_id(), plasmid.get_ice_id(), strain.get_ice_id()


def set_groups(filename, ice_url, ice_username, ice_password, group_name):
    '''Sets group permissions.'''
    with io.open(filename, 'r', newline='\r') as fle:
        ice_client = ice_utils.ICEClient(ice_url, ice_username, ice_password)
        groups = ice_client.get_groups()

        for ice_id in fle:
            ice_client.add_permission(ice_id, groups[group_name])


def main(args):
    '''main method.'''
    set_groups(*args)

if __name__ == '__main__':
    main(sys.argv[1:])
