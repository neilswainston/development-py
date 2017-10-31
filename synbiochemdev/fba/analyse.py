'''
development-py (c) University of Manchester 2017

development-py is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=too-few-public-methods
import sys

import cobra
from cobra.util.solver import linear_reaction_coefficients


class Analyser(object):
    '''Class to analyse a metabolic model.'''

    def __init__(self, filename):
        self.__model = cobra.io.read_sbml_model(filename)

    def analyse(self, obj='EX_pinocembrin', obj_low_bound_pc=0.8):
        '''Analyse model.'''

        solution = self.__model.optimize()
        print solution
        print

        # Fix objective lower bounds:
        initial_coeffs = linear_reaction_coefficients(self.__model)

        for coeff in initial_coeffs:
            coeff.lower_bound = solution.fluxes[coeff.id] * obj_low_bound_pc

        self.__model.objective = obj
        solution = self.__model.optimize()
        print solution
        print

        self.__model.summary()
        print

        for initial_coeff in initial_coeffs:
            print initial_coeff.id + '\t' + \
                str(solution.fluxes[initial_coeff.id])

        print

        fluxes = solution.fluxes[abs(solution.fluxes) >
                                 solution.objective_value / 100].to_frame()

        fluxes = fluxes.reset_index()
        fluxes['reaction'] = fluxes['index'].apply(lambda x:
                                                   self.__get_react_str(x))

        fluxes.sort_values(by='fluxes').to_csv('fluxes.csv', index=False)

    def __get_react_str(self, react_id):
        '''Get reaction string from reaction id.'''
        reaction = self.__model.reactions.get_by_id(react_id)
        return reaction.build_reaction_string()


def main(args):
    '''main method.'''
    analyser = Analyser(args[0])
    analyser.analyse()


if __name__ == '__main__':
    main(sys.argv[1:])
