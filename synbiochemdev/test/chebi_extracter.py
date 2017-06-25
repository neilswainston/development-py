'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import csv
import sys

from libchebipy._chebi_entity import ChebiEntity


def extract(chebi_id, filename):
    '''Recursively extracts data from ChEBI ontology.'''
    chebi_data = {}
    _get_chebi_data(chebi_id, chebi_data)
    _write(chebi_data, filename)


def _get_chebi_data(chebi_id, chebi_data):
    '''Recursively extracts data from supplied ChEBI id.'''
    if chebi_id not in chebi_data:
        chebi_entity = ChebiEntity(chebi_id)
        hmdb_ids = [db_acc.get_accession_number()
                    for db_acc in chebi_entity.get_database_accessions()
                    if db_acc.get_source() == 'HMDB']
        chebi_data[chebi_id] = {'Name': chebi_entity.get_name(),
                                'ChEBI id': chebi_entity.get_id(),
                                'HMDB id': ';'.join(hmdb_ids)}

        for incoming in chebi_entity.get_incomings():
            _get_chebi_data(incoming.get_target_chebi_id(), chebi_data)


def _write(chebi_data, filename):
    '''Writes ChEBI data to file.'''
    with open(filename, 'wb') as fle:
        # Create the headings:
        writer = csv.DictWriter(
            fle, chebi_data.values()[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(chebi_data.values())


def main(argv):
    '''main method'''
    extract(*argv)


if __name__ == '__main__':
    main(sys.argv[1:])
