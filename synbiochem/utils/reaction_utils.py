'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
from scipy.optimize import linprog

import synbiochem.utils.chemistry_utils as chem_utils


def balance(reaction_def, optional_formulae=None,
            optional_charges=None, max_stoich=8):
    '''Applies linear programming to balance reaction.'''
    if optional_formulae is None:
        optional_formulae = ['H', 'H2O']

    if optional_charges is None:
        optional_charges = [1, 0]

    # [['CO2', 'C5H7O4'], ['C3H3O3'], ['H', 'H2O'], ['H', 'H2O']]
    all_formulae = [[x[0] for x in reaction_def if x[2] <= 0],
                    [x[0] for x in reaction_def if x[2] > 0]]
    all_formulae.extend([optional_formulae for _ in range(2)])

    # [[0, -1], [-1], [1, 0], [1, 0]]
    all_charges = [[x[1] for x in reaction_def if x[2] <= 0],
                   [x[1] for x in reaction_def if x[2] > 0]]
    all_charges.extend([optional_charges for _ in range(2)])

    all_elem_comp = [[_get_elem_comp(formula, idx)
                      for formula in formulae]
                     for idx, formulae in enumerate(all_formulae)]

    (balanced, stoichs) = _optimise(_get_elem_matrix(all_elem_comp,
                                                     all_charges),
                                    max_stoich, len(optional_formulae))

    return balanced, _format_stoichs(stoichs, all_formulae, all_charges) \
        if balanced else None


def _get_elem_comp(formula, idx):
    '''Returns elemental compositions, multiplying stoichiometry by -1
    if reactants.'''
    elem_comp = chem_utils.get_elem_comp(formula)
    elem_comp.update((x, y * (-1 if idx % 2 == 0 else 1))
                     for x, y in elem_comp.items())
    return elem_comp


def _get_elem_matrix(all_elem_comp, all_charges):
    '''Gets the elemental (and charge) matrix representing rows of elements
    and columns of compounds.'''
    a_matrix = []
    elements = [elem_comp.keys() for elem_comps in all_elem_comp
                for elem_comp in elem_comps]

    for element in set(itertools.chain(*elements)):
        a_matrix.append([elem_comp[element] if element in elem_comp else 0
                         for elem_comps in all_elem_comp
                         for elem_comp in elem_comps])

    a_matrix.append(list(itertools.chain(*([charge * (-1 if idx % 2 == 0
                                                      else 1)
                                            for charge in charges]
                                           for idx, charges
                                           in enumerate(all_charges)))))

    return a_matrix


def _optimise(a_matrix, max_stoich, num_opt_formulae):
    '''Optimise linear program and return result.'''
    bounds = [(1, max_stoich)] * (len(a_matrix[0]) - num_opt_formulae * 2) + \
        [(0, max_stoich)] * num_opt_formulae * 2

    res = linprog([1] * len(a_matrix[0]),
                  A_eq=a_matrix,
                  b_eq=[0] * (len(a_matrix)),
                  bounds=bounds)

    return res.success, res.x


def _format_stoichs(stoichs, all_formulae, all_charges):
    '''Formats the input into (formula, charge, stoichiometry) tuples.'''
    return [(a[0], b, a[1] * c)
            for a, b, c in zip([(x, -1 if idx % 2 == 0 else 1)
                                for idx, formulae in enumerate(all_formulae)
                                for x in formulae],
                               [y for charges in all_charges for y in charges],
                               stoichs)
            if c > 1e-8]


def main():
    '''main method.'''
    print balance([('CO2', 0, -1.0), ('C5H7O4', -1, -1.0),
                   ('C3H3O3', -1, 1.0)])

if __name__ == '__main__':
    main()
