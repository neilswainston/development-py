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

    def __init__(self, post_seq, len_target=50, tir_target=None,
                 dg_target=None):
        # Check if dg_total or TIR (translation initiation rate) was specified.
        # If TIR, then convert to dg_total.
        self.__dg_target = _RBS_CALC.RT_eff * \
            (_RBS_CALC.logK - math.log(float(tir_target))) \
            if tir_target is not None else dg_target

        # If an initial RBS (sequences[1] is given, use it.
        # Otherwise, randomly choose one that is a decent starting point.
        (rbs, _) = RBS_MC_Design.GetInitialRBS('', post_seq,
                                               self.__dg_target)

        pre_seq = \
            ''.join([random.choice(['A', 'T', 'G', 'C'])
                     for _ in range(0, max(0, len_target - len(rbs)))])

        self.__sequences = [pre_seq, rbs, post_seq]
        self.__dg = self.__calc_dg(rbs)
        self.__sequences_new = [None, None, post_seq]
        self.__dg_new = None

    def get_energy(self):
        '''Gets the (simulated annealing) energy.'''
        return self.__dg - self.__dg_target

    def get_dg(self):
        '''Gets the delta G.'''
        return self.__dg

    def get_tir(self):
        '''Gets the translation initiation rate.'''
        return _RBS_CALC.K * math.exp(-self.__dg / _RBS_CALC.RT_eff)

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
            pre_seq_new = self.__sequences[0][1:] \
                if self.__sequences[0] > 0 else ''

        elif move == 'delete' and len(self.__sequences[1]) > 1:
            rbs_new = self.__sequences[1][0:pos] + \
                self.__sequences[1][pos + 1:len(self.__sequences[1])]
            pre_seq_new = random.choice(['A', 'T', 'G', 'C']) + \
                self.__sequences[0]

        elif move == 'replace':
            letter = random.choice(['A', 'T', 'G', 'C'])
            rbs_new = self.__sequences[1][0:pos] + letter + \
                self.__sequences[1][pos + 1:len(self.__sequences[1])]
            pre_seq_new = self.__sequences[0]

        else:
            pre_seq_new = self.__sequences[0]
            rbs_new = self.__sequences[1]

        rbs_new = RBS_MC_Design.RemoveStartCodons(rbs_new)

        self.__sequences_new[0] = pre_seq_new
        self.__sequences_new[1] = rbs_new
        self.__dg_new = self.__calc_dg(rbs_new, verbose)
        return self.__dg_new - self.__dg_target

    def accept(self):
        '''Accept potential update.'''
        self.__sequences = self.__sequences_new
        self.__dg = self.__dg_new
        self.__sequences_new = [None, None, self.__sequences[2]]
        self.__dg_new = None

    def __calc_dg(self, rbs, verbose=False):
        '''Calculates (simulated annealing) energy for given RBS.'''
        calc = RBS_MC_Design.Run_RBS_Calculator(self.__sequences[0],
                                                self.__sequences[2],
                                                rbs,
                                                verbose)

        return calc.dG_total_list[0]

    def __repr__(self):
        # return '%r' % (self.__dict__)
        return str(self.__dg) + '\t' + str(self.get_tir()) + '\t' + \
            self.__sequences[0] + ' ' + \
            self.__sequences[1] + ' ' + self.__sequences[2]

    def __print__(self):
        return self.__repr__


def main(argv):
    '''main method.'''
    print sim_ann.optimise(RBSSolution(argv[1], len_target=int(argv[2]),
                                       tir_target=float(argv[3])),
                           verbose=False)


if __name__ == '__main__':
    main(sys.argv)
