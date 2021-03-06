'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import json
import sys
import urllib2

import synbiochem.utils.chem_utils as chem_utils


def get_metabolome(model_ids):
    '''Gets the metabolome from supplied BiGG model ids.'''
    bigg_data = {}
    organisms = {}
    databases = set()

    for model_id in model_ids:
        model_url = 'http://bigg.ucsd.edu/api/v2/models/' + model_id + \
            '/metabolites'
        model_data = json.loads(urllib2.urlopen(model_url).read())

        for current_bigg_ids in [(result['bigg_id'],
                                  result['bigg_id'] + '_' +
                                  result['compartment_bigg_id'],
                                  result['organism'])
                                 for result in model_data['results']]:

            organisms[model_id] = current_bigg_ids[2]

            if current_bigg_ids[0] not in bigg_data:
                databases.update(_add_metabolite(model_id, model_ids,
                                                 current_bigg_ids,
                                                 bigg_data))

    return bigg_data, organisms, databases


def _add_metabolite(model_id, model_ids, current_bigg_ids,
                    bigg_data):
    '''Adds a given metabolite to the collection of metabolites.'''
    met_url = 'http://bigg.ucsd.edu/api/v2/models/' + model_id + \
        '/metabolites/' + current_bigg_ids[1]
    met_data = json.loads(urllib2.urlopen(met_url).read())
    bg_dat = {}
    bigg_data[met_data['bigg_id']] = bg_dat
    bg_dat['name'] = met_data['name']
    bg_dat['formula'] = met_data['formula']
    bg_dat['database_links'] = dict((key, value[0]['id'])
                                    for key, value
                                    in met_data['database_links']
                                    .iteritems())

    bg_dat['models'] = [model['bigg_id']
                        for model
                        in met_data['other_models_with_metabolite']
                        if model['bigg_id'] in model_ids]
    bg_dat['models'].append(model_id)
    bg_dat['models'] = sorted(set(bg_dat['models']))

    try:
        bg_dat['mass'] = \
            chem_utils.get_molecular_mass(met_data['formula'])
    except KeyError, error:
        print met_data['bigg_id'] + ': ' + str(error)

    return bg_dat['database_links'].keys()


def main(argv):
    '''main method.'''
    bigg_data, organisms, databases = get_metabolome(argv[1:])
    chem_props = ['name', 'formula', 'mass']
    headers = ['id'] + chem_props + list(databases)

    print '\t'.join(headers + organisms.keys())
    print '\t'.join(['' for _ in range(len(headers))] + organisms.values())

    for bigg_id in bigg_data:
        print '\t'.join([bigg_id] + [str(bigg_data[bigg_id].get(chem_prop))
                                     for chem_prop in chem_props] +
                        [bigg_data[bigg_id]['database_links'][database]
                         if database in bigg_data[bigg_id]['database_links']
                         else ''
                         for database in list(databases)] +
                        [str(organism in bigg_data[bigg_id].get('models'))
                         for organism in organisms])


if __name__ == '__main__':
    main(sys.argv)
