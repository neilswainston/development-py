'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochem.build.rbs import RBSSolution
import synbiochem.optimisation.simulated_annealing as sim_ann


def main():
    '''main method.'''
    taxonomy_id = '83333'
    len_target = 60
    tir_target = 15000

    optimised_chi = 'MGSMALPSISAKWKGKNAKELTESVPFFRQLVTGEFEKLARVTMKKRLTGIQY' + \
        'SEKVVENCEEIMKASGKYTRSEAKAIDQFLMVFKNQDFPPGSSIIFAICPKGSLTIAFSKEERVP' + \
        'KTGKAVIKNKLLGEAVLESMIGKNGVSPATRKSLAERLSKLMNKKDPYNEANVNVATKN'

    # TAL
    sol = sim_ann.optimise(RBSSolution(['Q3IWB0', 'A9AUJ9', 'A5FKY3', 'Q1LRV9',
                                        'Q8VXG7'], taxonomy_id, len_target,
                                       tir_target), verbose=True)[0]
    sol.print_sol()

    # PAL
    sol = sim_ann.optimise(RBSSolution(['P35510'], taxonomy_id, len_target,
                                       tir_target))[0]
    sol.print_sol()

    # 4CL
    sol = sim_ann.optimise(RBSSolution(['O54075', 'O24146', 'Q9S777', 'P31687',
                                        'Q42879', 'Q9K3W1'], taxonomy_id,
                                       len_target, tir_target),
                           verbose=True)[0]
    sol.print_sol()

    # CHS
    sol = sim_ann.optimise(RBSSolution(['P13114'], taxonomy_id, len_target,
                                       tir_target))[0]
    sol.print_sol()

    # CHI
    sol = sim_ann.optimise(RBSSolution(['P41088', optimised_chi], taxonomy_id,
                                       len_target, tir_target))[0]
    sol.print_sol()

if __name__ == '__main__':
    main()
