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
import synbiochem.utils.sequence_utils as seq_util
import synbiochem.optimisation.simulated_annealing as sim_ann


# Necessary to get constants hidden as class variables in RBS_Calculator:
_RBS_CALC = RBS_Calculator.RBS_Calculator('A', [0, 0])


class RBSSolution(object):
    '''Solution for RBS optimisation.'''

    def __init__(self, prot_seq, taxonomy_id, len_target=50, tir_target=None):
        # Check if dg_total or TIR (translation initiation rate) was specified.
        # If TIR, then convert to dg_total.
        self.__dg_target = _RBS_CALC.RT_eff * \
            (_RBS_CALC.logK - math.log(float(tir_target)))

        self.__prot_seq = prot_seq
        self.__cod_opt = seq_util.CodonOptimiser(taxonomy_id)
        cds = self.__cod_opt.get_codon_optimised_seq(self.__prot_seq)

        # If an initial RBS (sequences[1] is given, use it.
        # Otherwise, randomly choose one that is a decent starting point.
        (rbs, _) = RBS_MC_Design.GetInitialRBS('', cds,
                                               self.__dg_target)

        pre_seq = \
            ''.join([random.choice(['A', 'T', 'G', 'C'])
                     for _ in range(0, max(0, len_target - len(rbs)))])

        self.__sequences = [pre_seq, rbs, cds]
        self.__dg = self.__calc_dg(rbs)
        self.__sequences_new = [None, None, cds]
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
        '''Mutates and scores whole design.'''
        self.__mutate_pre_seq()
        self.__mutate_rbs()
        self.__mutate_cds()
        self.__dg_new = self.__calc_dg(self.__sequences_new[1], verbose)
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

    def __mutate_pre_seq(self):
        '''Mutates pre-sequence.'''
        pos = int(random.random() * len(self.__sequences[0]))
        self.__sequences_new[0] = _replace(self.__sequences[0], pos,
                                           _rand_nuc())

    def __mutate_rbs(self):
        '''Mutates RBS.'''
        weighted_moves = [('insert', 0.1), ('delete', 0.1), ('replace', 0.8)]
        move = RBS_MC_Design.weighted_choice(weighted_moves)
        pos = int(random.random() * len(self.__sequences[1]))

        if move == 'insert' and \
                len(self.__sequences[1]) < RBS_MC_Design.Max_RBS_Length:
            letter = random.choice(['A', 'T', 'G', 'C'])
            rbs_new = self.__sequences[1][0:pos] + letter + \
                self.__sequences[1][pos:len(self.__sequences[1])]
            pre_seq_new = self.__sequences_new[0][1:] \
                if len(self.__sequences_new[0]) > 0 else ''

        elif move == 'delete' and len(self.__sequences[1]) > 1:
            rbs_new = _replace(self.__sequences[1], pos, '')
            pre_seq_new = random.choice(['A', 'T', 'G', 'C']) + \
                self.__sequences_new[0]

        elif move == 'replace':
            rbs_new = _replace(self.__sequences[1], pos, _rand_nuc())
            pre_seq_new = self.__sequences_new[0]

        else:
            pre_seq_new = self.__sequences_new[0]
            rbs_new = self.__sequences[1]

        rbs_new = RBS_MC_Design.RemoveStartCodons(rbs_new)

        self.__sequences_new[0] = pre_seq_new
        self.__sequences_new[1] = rbs_new

    def __mutate_cds(self):
        '''Mutates CDS.'''
        self.__sequences_new[2] = \
            self.__cod_opt.mutate(self.__prot_seq, self.__sequences[2],
                                  1.0/len(self.__sequences[2]))

    def __repr__(self):
        # return '%r' % (self.__dict__)
        return str(self.__dg) + '\t' + str(self.get_tir()) + '\t' + \
            self.__sequences[0] + ' ' + \
            self.__sequences[1] + ' ' + self.__sequences[2]

    def __print__(self):
        return self.__repr__


def _replace(sequence, pos, nuc):
    '''Replace nucleotide at pos with nuc.'''
    return sequence[:pos] + nuc + sequence[pos + 1:]


def _rand_nuc():
    '''Returns a random nucleotide.'''
    return random.choice(['A', 'T', 'G', 'C'])


def main(argv):
    '''main method.'''
    print sim_ann.optimise(RBSSolution(argv[1], argv[2],
                                       len_target=int(argv[3]),
                                       tir_target=float(argv[4])),
                           verbose=False)


if __name__ == '__main__':
    main(sys.argv)
