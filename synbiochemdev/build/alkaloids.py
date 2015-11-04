'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochemdev.build.rbs import RBSSolution
import synbiochemdev.optimisation.simulated_annealing as sim_ann


def main():
    '''main method.'''
    taxonomy_id = '83333'
    len_target = 60
    tir_target = 15000
    verbose = True

    optimised_p0ab91 = 'MNYQNDDLRIKEIKELLPPVALLEKFPATENAANTVAHARKAIHKILKGN' + \
        'DDRLLVVIGPCSIHDPVAAKEYATRLLALREELKDELEIVMRVYFEKPRTTVGWKGLINDPHMDN' + \
        'SFQINDGLRIARKLLLDINDSGLPAAGEFLNMITPQYLADLMSWGAIGARTTESQVHRELASGLS' + \
        'CPVGFKNGTDGTIKVAIDAINAAGAPHCFLSVTKWGHSAIVNTSGNGDCHIILRGGKEPNYSAKH' + \
        'VAEVKEGLNKAGLPAQVMIDFSHANSSKQFKKQMDVCADVCQQIAGGEKAIIGVMVESHLVEGNQ' + \
        'SLESGEPLAYGKSITDACIGWEDTDALLRQLANAVKARRG'

    sol = sim_ann.optimise(RBSSolution(['P07023'], taxonomy_id, len_target,
                                       tir_target), verbose=verbose)[0]
    sol.print_sol()

    sol = sim_ann.optimise(RBSSolution([optimised_p0ab91], taxonomy_id,
                                       len_target, tir_target),
                           verbose=verbose)[0]
    sol.print_sol()

    sol = sim_ann.optimise(RBSSolution(['P23538'], taxonomy_id, len_target,
                                       tir_target), verbose=verbose)[0]
    sol.print_sol()

    sol = sim_ann.optimise(RBSSolution(['P27302'], taxonomy_id, len_target,
                                       tir_target), verbose=verbose)[0]
    sol.print_sol()

    sol = sim_ann.optimise(RBSSolution(['Q57160'], taxonomy_id, len_target,
                                       tir_target), verbose=verbose)[0]
    sol.print_sol()

    sol = sim_ann.optimise(RBSSolution(['Q57501'], taxonomy_id, len_target,
                                       tir_target), verbose=verbose)[0]
    sol.print_sol()

    sol = sim_ann.optimise(RBSSolution(['K0IQX2'], taxonomy_id, len_target,
                                       tir_target), verbose=verbose)[0]
    sol.print_sol()

if __name__ == '__main__':
    main()
