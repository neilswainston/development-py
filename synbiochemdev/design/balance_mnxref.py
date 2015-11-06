'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochem.design.mnxref import MnxRefReader
import synbiochem.utils.chem_utils as chem_utils


def main():
    '''Main method.'''
    reader = MnxRefReader()
    chem_data = reader.get_chem_data()

    for reac_id, reac_props in reader.get_reac_data().iteritems():
        participants = chem_utils.parse_equation(reac_props['equation'])

        reac_def = [(chem_data[prt[0]]['formula'], chem_data[prt[0]]['charge'],
                     prt[1])
                    for prt in participants
                    if 'formula' in chem_data[prt[0]] and 'charge'
                    in chem_data[prt[0]]]

        is_balanced, was_balanced, balanced_def = chem_utils.balance(reac_def)

        print '\t'.join([reac_id, str(is_balanced), str(was_balanced),
                         reac_props['description'], str(balanced_def)])

if __name__ == '__main__':
    main()
