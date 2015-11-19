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

    for chem_id, props in chem_data.iteritems():
        if 'formula' not in props:
            print chem_id + '\t' + str(props)

    for reac_id, reac_props in reader.get_reac_data().iteritems():
        participants = chem_utils.parse_equation(reac_props['equation'])

        try:
            reac_def = [(chem_data[prt[0]]['formula'],
                         chem_data[prt[0]]['charge']
                         if 'charge' in chem_data[prt[0]] else 0,
                         prt[1])
                        for prt in participants]

            is_balanced, was_balanced, balanced_def = \
                chem_utils.balance(reac_def)

            print '\t'.join([reac_id, str(is_balanced),
                             str(was_balanced),
                             reac_props['description'],
                             '' if balanced_def is None
                             else str(balanced_def)])
        except KeyError:
            print '\t'.join([reac_id, 'UNKNOWN', 'UNKNOWN',
                             reac_props['description'], ''])

if __name__ == '__main__':
    main()
