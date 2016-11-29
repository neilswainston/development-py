'''
Created on 29 Nov 2016

@author: neilswainston
'''
import io
import sys

from synbiochem.utils import ice_utils


def main(args):
    '''main method.'''
    with io.open(args[0], 'r', newline='\r') as fle:
        pairs = [str(line).split() for line in fle]

        ice_client = ice_utils.ICEClient(args[1], args[2], args[3])

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

            print host.get_ice_id(), plasmid.get_ice_id(), strain.get_ice_id()


if __name__ == '__main__':
    main(sys.argv[1:])
