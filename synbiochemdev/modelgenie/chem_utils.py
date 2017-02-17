'''
ModelGenie (c) University of Manchester 2015

ModelGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import re


def parse_equation(equation, separator):
    '''Parses a chemical equation, returning the participants
    as a list of name, stoichiometry tuples.
    Negative stoichiometries define reactants, positive products.'''
    equation_terms = [re.split('\\s*\\+\\s*', equation_side)
                      for equation_side in
                      re.split('\\s*' + separator + '\\s*', equation)]

    participants = []

    # Add reactants:
    __add_reaction_participants(equation_terms[0], -1, participants)

    # Add products:
    __add_reaction_participants(equation_terms[1], 1, participants)

    return participants


def __add_reaction_participants(equation_term, stoich_factor, participants):
    '''Adds participants to reaction.'''
    for participant in equation_term:
        terms = participant.split()
        participants.append((participant, stoich_factor) if len(terms) == 1
                            else (terms[1], stoich_factor * float(terms[0])))
