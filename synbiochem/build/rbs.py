'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math
import random
import sys

import RBS_Calculator
import RBS_MC_Design
import synbiochem.optimisation.simulated_annealing as sim_ann


# Necessary to get constants hidden as class variables in RBS_Calculator:
_RBS_CALC = RBS_Calculator.RBS_Calculator('A', [0, 0])


class RBSSolution(object):
    '''Solution for RBS optimisation.'''

    def __init__(self, sequences, tir_target=None,
                 dg_target=None):
        # Check if dg_total or TIR (translation initiation rate) was specified.
        # If TIR, then convert to dg_total.
        self.__dg_target = _RBS_CALC.RT_eff * \
            (_RBS_CALC.logK - math.log(float(tir_target))) \
            if tir_target is not None else dg_target

        self.__sequences = sequences

        # If an initial RBS (sequences[1] is given, use it.
        # Otherwise, randomly choose one that is a decent starting point.
        if sequences[1] is None:
            (sequences[1], _) = RBS_MC_Design.GetInitialRBS(sequences[0],
                                                            sequences[2],
                                                            self.__dg_target)

        self.__energy = self.__calc_energy(sequences[1])
        self.__rbs_new = None
        self.__energy_new = None

    def get_energy(self):
        '''Gets the (simulated annealing) energy.'''
        return self.__energy

    def mutate(self, verbose=False):
        '''Mutates and scores RBS.'''
        weighted_moves = [('insert', 0.1), ('delete', 0.1), ('replace', 0.8)]
        move = RBS_MC_Design.weighted_choice(weighted_moves)
        pos = int(random.random() * len(self.__sequences[1]))

        if move == 'insert' and \
                len(self.__sequences[1]) < RBS_MC_Design.Max_RBS_Length:
            letter = random.choice(['A', 'T', 'G', 'C'])
            rbs_new = self.__sequences[1][0:pos] + letter + \
                self.__sequences[1][pos:len(self.__sequences[1])]
        elif move == 'delete' and len(self.__sequences[1]) > 1:
            rbs_new = self.__sequences[1][0:pos] + \
                self.__sequences[1][pos + 1:len(self.__sequences[1])]
        elif move == 'replace':
            letter = random.choice(['A', 'T', 'G', 'C'])
            rbs_new = self.__sequences[1][0:pos] + letter + \
                self.__sequences[1][pos + 1:len(self.__sequences[1])]
        else:
            rbs_new = self.__sequences[1]

        self.__rbs_new = RBS_MC_Design.RemoveStartCodons(rbs_new)
        self.__energy_new = self.__calc_energy(rbs_new, verbose)
        return self.__energy_new

    def accept(self):
        '''Accept potential update.'''
        self.__sequences[1] = self.__rbs_new
        self.__energy = self.__energy_new
        self.__rbs_new = None
        self.__energy_new = None

    def __calc_energy(self, rbs, verbose=False):
        '''Calculates (simulated annealing) energy for given RBS.'''
        calc = RBS_MC_Design.Run_RBS_Calculator(self.__sequences[0],
                                                self.__sequences[2],
                                                rbs,
                                                verbose)

        return calc.dG_total_list[0] - self.__dg_target

    def __repr__(self):
        # return '%r' % (self.__dict__)
        return self.__sequences[1]

    def __print__(self):
        return self.__repr__


def main(argv):
    '''main method.'''
    print sim_ann.optimise(RBSSolution([argv[1], None, argv[2]],
                                       tir_target=float(argv[3])),
                           verbose=False)


if __name__ == '__main__':
    main(sys.argv)
