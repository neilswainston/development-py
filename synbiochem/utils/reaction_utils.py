'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import itertools
from scipy.optimize import linprog

import chemistry_utils


def balance():
    '''balance method.'''
    ccc = [-1, 4]
    aaa = [[-3, 1], [1, 2]]
    bbb = [6, 4]
    x0_bounds = (None, None)
    x1_bounds = (-3, None)
    res = linprog(ccc, A_ub=aaa, b_ub=bbb,
                  bounds=(x0_bounds, x1_bounds),
                  options={'disp': True})
    print res


def balance2():
    '''balance method.'''
    ccc = [1, 1, 1, 1, 1, 1, 1, 1]
    aaa = [[1, 5, -3, 0, 0, 0, 0, -1],
           [2, 4, -3, 0, 1, 0, -1, -2],
           [0, 7, -3, 1, 2, -1, -2, 0],
           [0, -1, 1, 1, 0, -1, 0, 0]]

    bbb = [0, 0, 0, 0]
    bounds = ((1, 8), (1, 8), (1, 8), (0, 8), (0, 8), (0, 8), (0, 8), (0, 8))
    res = linprog(ccc, A_ub=aaa, b_ub=bbb, bounds=bounds,
                  options={'disp': True})
    print res


def balance3(all_formulae, all_charges):
    optional_formulae = ['H', 'H2O']
    optional_charges = [1, 0]

    all_formulae.extend([optional_formulae for _ in range(2)])
    all_charges.extend([optional_charges for _ in range(2)])
    all_elem_comp = [[_get_elem_comp(formula, idx)
                      for formula in formulae]
                     for idx, formulae in enumerate(all_formulae)]

    elements = [elem_comp.keys() for elem_comps in all_elem_comp
                for elem_comp in elem_comps]

    aaa = []
    for element in set(itertools.chain(*elements)):
        aaa.append([elem_comp[element] if element in elem_comp else 0
                    for elem_comps in all_elem_comp
                    for elem_comp in elem_comps])

    aaa.append(list(itertools.chain(*([charge * (-1 if idx % 2 == 0 else 1)
                                       for charge in charges]
                                      for idx,
                                      charges in enumerate(all_charges)))))

    ccc = [1] * len(elements)
    bbb = [0] * len(aaa)
    bounds = ((1, 4), (1, 4), (1, 4), (0, 4), (0, 4), (0, 4), (0, 4))
    res = linprog(ccc, A_ub=aaa, b_ub=bbb, bounds=bounds,
                  options={'disp': True})
    print res

def _get_elem_comp(formula, idx):
    '''Returns elemental compositions, multiplying stoichiometry by -1
    if reactants.'''
    elem_comp = chemistry_utils.get_elem_comp(formula)
    elem_comp.update((x, y * (-1 if idx % 2 == 0 else 1))
                     for x, y in elem_comp.items())
    return elem_comp


def main():
    '''main method.'''
    balance3([['CO2', 'C5H7O4'], ['C3H3O3']], [[0, -1], [-1]])

if __name__ == '__main__':
    main()
