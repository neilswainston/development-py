'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochem.build.rbs import RBSSolution
import synbiochem.optimisation.simulated_annealing as sim_ann


if __name__ == '__main__':
    '''main method.'''
    taxonomy_id = '83333'
    len_target = 60
    tir_target = 15000

    optimised_P0AB91 = 'MNYQNDDLRIKEIKELLPPVALLEKFPATENAANTVAHARKAIHKILKGN' + \
        'DDRLLVVIGPCSIHDPVAAKEYATRLLALREELKDELEIVMRVYFEKPRTTVGWKGLINDPHMDN' + \
        'SFQINDGLRIARKLLLDINDSGLPAAGEFLNMITPQYLADLMSWGAIGARTTESQVHRELASGLS' + \
        'CPVGFKNGTDGTIKVAIDAINAAGAPHCFLSVTKWGHSAIVNTSGNGDCHIILRGGKEPNYSAKH' + \
        'VAEVKEGLNKAGLPAQVMIDFSHANSSKQFKKQMDVCADVCQQIAGGEKAIIGVMVESHLVEGNQ' + \
        'SLESGEPLAYGKSITDACIGWEDTDALLRQLANAVKARRG'

    print sim_ann.optimise(RBSSolution(['P07023'], taxonomy_id, len_target,
                                       tir_target), verbose=True)

    print sim_ann.optimise(RBSSolution([optimised_P0AB91], taxonomy_id,
                                       len_target, tir_target))

    print sim_ann.optimise(RBSSolution(['P23538'], taxonomy_id, len_target,
                                       tir_target))

    print sim_ann.optimise(RBSSolution(['P27302'], taxonomy_id, len_target,
                                       tir_target))

    print sim_ann.optimise(RBSSolution(['Q57160'], taxonomy_id, len_target,
                                       tir_target))

    print sim_ann.optimise(RBSSolution(['Q57501'], taxonomy_id, len_target,
                                       tir_target))

    print sim_ann.optimise(RBSSolution(['K0IQX2'], taxonomy_id, len_target,
                                       tir_target))
