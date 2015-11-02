'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods

import os
import re


class MolecularMassCalculator(object):
    '''Class for calculating molecular masses.'''

    def __init__(self):
        self.__element_to_mass = {}
        location = \
            os.path.realpath(os.path.join(os.getcwd(),
                                          os.path.dirname(__file__)))

        with open(os.path.join(location, 'element_data.txt')) as input_file:
            for line in input_file.readlines():
                (element, mass) = line.strip().split()
                self.__element_to_mass[element] = float(mass)

    def get_molecular_mass(self, formula):
        '''Calculate and return molecular mass from chemical formula.'''
        return sum([self.__element_to_mass[element] * count
                    for element, count in get_elem_comp(formula)])


def get_elem_comp(formula):
    '''Gets elemental composition as a dict from formula.'''
    elem_comp = {}

    for term in re.findall('[A-Z]{1}[0-9]*[a-z]{0,1}[0-9]*', formula):
        element = re.search('[A-z]*', term).group(0)
        result = re.search('[0-9]+', term)
        elem_comp[element] = int(result.group(0)) if result else 1

    return elem_comp
