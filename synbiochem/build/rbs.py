'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math
import random
import re
import sys

import RBS_Calculator
import RBS_MC_Design
import synbiochem.utils.sequence_utils as seq_utils
import synbiochem.utils.uniprot_utils as uniprot_utils
import synbiochem.optimisation.simulated_annealing as sim_ann


# Necessary to get constants hidden as class variables in RBS_Calculator:
_RBS_CALC = RBS_Calculator.RBS_Calculator('A', [0, 0])


class RBSSolution(object):
    '''Solution for RBS optimisation.'''

    def __init__(self, uniprot_id_set, taxonomy_id, len_target, tir_target):
        # Check if dg_total or TIR (translation initiation rate) was specified.
        # If TIR, then convert to dg_total.
        self.__dg_target = _RBS_CALC.RT_eff * \
            (_RBS_CALC.logK - math.log(float(tir_target)))

        self.__prot_seqs = uniprot_utils.get_sequences(uniprot_id_set)
        self.__cod_opt = seq_utils.CodonOptimiser(taxonomy_id)
        cds = [self.__cod_opt.get_codon_optimised_seq(prot_seq)
               for prot_seq in self.__prot_seqs.values()]

        # Randomly choose an RBS that is a decent starting point,
        # using the first CDS as the upstream sequence:
        (rbs, _) = RBS_MC_Design.GetInitialRBS('', cds[0],
                                               self.__dg_target)

        pre_seq = \
            ''.join([random.choice(['A', 'T', 'G', 'C'])
                     for _ in range(0, max(0, len_target - len(rbs)))])

        self.__seqs = [pre_seq, rbs, cds]
        self.__dgs = self.__calc_dgs(rbs)
        self.__seqs_new = [None, None, cds]
        self.__dgs_new = None

    def get_energy(self, dgs=None, cdss=None):
        '''Gets the (simulated annealing) energy.'''
        dgs = self.__dgs if dgs is None else dgs
        cdss = self.__seqs[2] if cdss is None else cdss
        cais = [self.__cod_opt.get_cai(cds) for cds in cdss]
        return sum([abs(d_g-self.__dg_target) for d_g in dgs])/len(dgs) * \
            (1-sum(cais)/len(cais))

    def mutate(self, verbose=False):
        '''Mutates and scores whole design.'''
        self.__mutate_pre_seq()
        self.__mutate_rbs()
        self.__mutate_cds()
        self.__dgs_new = self.__calc_dgs(self.__seqs_new[1], verbose)
        return self.get_energy(self.__dgs_new, self.__seqs_new[2])

    def accept(self):
        '''Accept potential update.'''
        self.__seqs = self.__seqs_new
        self.__dgs = self.__dgs_new
        self.__seqs_new = [None, None, self.__seqs[2]]
        self.__dgs_new = None

    def __calc_dgs(self, rbs, verbose=False):
        '''Calculates (simulated annealing) energies for given RBS.'''
        return [RBS_MC_Design.Run_RBS_Calculator(self.__seqs[0],
                                                 post_seq,
                                                 rbs,
                                                 verbose).dG_total_list[0]
                for post_seq in self.__seqs[2]]

    def __mutate_pre_seq(self):
        '''Mutates pre-sequence.'''
        pos = int(random.random() * len(self.__seqs[0]))
        self.__seqs_new[0] = _replace(self.__seqs[0], pos, _rand_nuc())

    def __mutate_rbs(self):
        '''Mutates RBS.'''
        weighted_moves = [('insert', 0.1), ('delete', 0.1), ('replace', 0.8)]
        move = RBS_MC_Design.weighted_choice(weighted_moves)
        pos = int(random.random() * len(self.__seqs[1]))

        if move == 'insert' and \
                len(self.__seqs[1]) < RBS_MC_Design.Max_RBS_Length:
            letter = random.choice(['A', 'T', 'G', 'C'])
            rbs_new = self.__seqs[1][0:pos] + letter + \
                self.__seqs[1][pos:len(self.__seqs[1])]
            pre_seq_new = self.__seqs_new[0][1:] \
                if len(self.__seqs_new[0]) > 0 else ''

        elif move == 'delete' and len(self.__seqs[1]) > 1:
            rbs_new = _replace(self.__seqs[1], pos, '')
            pre_seq_new = random.choice(['A', 'T', 'G', 'C']) + \
                self.__seqs_new[0]

        elif move == 'replace':
            rbs_new = _replace(self.__seqs[1], pos, _rand_nuc())
            pre_seq_new = self.__seqs_new[0]

        else:
            pre_seq_new = self.__seqs_new[0]
            rbs_new = self.__seqs[1]

        rbs_new = RBS_MC_Design.RemoveStartCodons(rbs_new)

        self.__seqs_new[0] = pre_seq_new
        self.__seqs_new[1] = rbs_new

    def __mutate_cds(self):
        '''Mutates CDS.'''
        self.__seqs_new[2] = \
            [self.__cod_opt.mutate(prot_seq, dna_seq, 1.0/len(dna_seq))
             for dna_seq, prot_seq in zip(self.__seqs[2],
                                          self.__prot_seqs.values())]

    def __repr__(self):
        # return '%r' % (self.__dict__)
        return str([self.__cod_opt.get_cai(prot_seq)
                    for prot_seq in self.__seqs[2]]) + '\t' + \
            str([_count_undesired_patterns(seq)
                 for seq in self.__seqs]) + '\t' + \
            str(_get_tirs(self.__dgs)) + '\t' + self.__seqs[0] + ' ' + \
            self.__seqs[1] + ' ' + str(self.__seqs[2])

    def __print__(self):
        return self.__repr__


def _get_tirs(dgs):
    '''Gets the translation initiation rate.'''
    return [_RBS_CALC.K * math.exp(-d_g / _RBS_CALC.RT_eff)
            for d_g in dgs]


def _replace(sequence, pos, nuc):
    '''Replace nucleotide at pos with nuc.'''
    return sequence[:pos] + nuc + sequence[pos + 1:]


def _rand_nuc():
    '''Returns a random nucleotide.'''
    return random.choice(['A', 'T', 'G', 'C'])


def _count_undesired_patterns(seqs):
    '''Counts undesired patterns in sequence.'''
    max_repeat_nucs = 4

    if isinstance(seqs, str) or isinstance(seqs, unicode):
        # Start codons | restriction sites | repeating nucleotides
        patterns = '|'.join(['[AGT]TG', 'GGTCTC', 'CACCTGC'] +
                            [x*max_repeat_nucs for x in ['A', 'C', 'G', 'T']])
        return len(re.findall(patterns, seqs))
    else:
        return [_count_undesired_patterns(seq) for seq in seqs]


def main(argv):
    '''main method.'''
    print sim_ann.optimise(RBSSolution(argv[4:], argv[1],
                                       len_target=int(argv[2]),
                                       tir_target=float(argv[3])),
                           verbose=False)


if __name__ == '__main__':
    main(sys.argv)
