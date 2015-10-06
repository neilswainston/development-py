'''
Created on 6 Oct 2015

@author: neilswainston
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
        mass = 0

        for term in re.findall('[A-Z]{1}[0-9]*[a-z]{0,1}[0-9]*',
                               formula):

            element = re.search("[A-z]*", term).group(0)
            result = re.search("[0-9]+", term)
            count = int(result.group(0)) if result else 1
            mass += self.__element_to_mass[element] * count

        return mass
