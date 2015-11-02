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
import synbiochem.utils.sequence_utils as seq_utils


# Necessary to get constants hidden as class variables in RBS_Calculator:
_RBS_CALC = RBS_Calculator.RBS_Calculator('A', [0, 0])

# Invalid pattern is restriction sites | repeating nucleotides:
_MAX_REPEATING_NUCS = 6
_INVALID_PATTERN = '|'.join(['GGTCTC', 'CACCTGC', 'GAATTC', 'AGATCT', 'GGATCC',
                             'CTCGAG', 'GAGTC[ACGT]{5}'] +
                            [x * _MAX_REPEATING_NUCS
                             for x in ['A', 'C', 'G', 'T']])


class RBSSolution(object):
    '''Solution for RBS optimisation.'''

    def __init__(self, protein_ids, taxonomy_id, len_target, tir_target):
        # Check if dg_total or TIR (translation initiation rate) was specified.
        # If TIR, then convert to dg_total.
        self.__dg_target = _RBS_CALC.RT_eff * \
            (_RBS_CALC.logK - math.log(float(tir_target)))

        self.__cod_opt = seq_utils.CodonOptimiser(taxonomy_id)

        self.__prot_seqs = seq_utils.get_sequences(protein_ids)
        cds = [self.__cod_opt.get_codon_optimised_seq(prot_seq,
                                                      _INVALID_PATTERN)
               for prot_seq in self.__prot_seqs.values()]

        stop_codon = self.__cod_opt.get_codon_optimised_seq('*')

        # Randomly choose an RBS that is a decent starting point,
        # using the first CDS as the upstream sequence:
        (rbs, _) = RBS_MC_Design.GetInitialRBS('', cds[0],
                                               self.__dg_target)

        post_seq_length = 30
        self.__seqs = [self.__get_valid_rand_seq(max(0,
                                                     len_target - len(rbs))),
                       rbs,
                       cds,
                       stop_codon,
                       self.__get_valid_rand_seq(post_seq_length)
                       if len(self.__prot_seqs) > 1 else None]
        self.__dgs = self.__calc_dgs(rbs)
        self.__seqs_new = [None, None, cds, self.__seqs[3], self.__seqs[4]]
        self.__dgs_new = None

    def get_energy(self, dgs=None, cdss=None):
        '''Gets the (simulated annealing) energy.'''
        dgs = self.__dgs if dgs is None else dgs
        cdss = self.__seqs[2] if cdss is None else cdss
        return sum([abs(d_g - self.__dg_target) for d_g in dgs]) / \
            len(self.__seqs[2]) * \
            (1 + (sum(seq_utils.count_pattern([self.__seqs[1]] + cdss,
                                              _INVALID_PATTERN)))**10)

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
        self.__seqs_new = [None, None, self.__seqs[2], self.__seqs[3],
                           self.__seqs[4]]
        self.__dgs_new = None

    def __calc_dgs(self, rbs, verbose=False):
        '''Calculates (simulated annealing) energies for given RBS.'''
        return [RBS_MC_Design.Run_RBS_Calculator(self.__seqs[0],
                                                 cds,
                                                 rbs,
                                                 verbose).dG_total_list[0]
                for cds in self.__seqs[2]]

    def __mutate_pre_seq(self):
        '''Mutates pre-sequence.'''
        pos = int(random.random() * len(self.__seqs[0]))
        pre_seq_new = _replace(self.__seqs[0], pos, _rand_nuc())

        if seq_utils.count_pattern(pre_seq_new + self.__seqs[1],
                                   _INVALID_PATTERN) + \
                seq_utils.count_pattern(pre_seq_new + self.__seqs[1],
                                        '[AGT]TG') == 0:
            self.__seqs_new[0] = pre_seq_new
        else:
            self.__seqs_new[0] = self.__seqs[0]

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

        if seq_utils.count_pattern(pre_seq_new + rbs_new, _INVALID_PATTERN) + \
                seq_utils.count_pattern(pre_seq_new + rbs_new, '[AGT]TG') == 0:
            self.__seqs_new[0] = pre_seq_new
            self.__seqs_new[1] = rbs_new
        else:
            self.__seqs_new[0] = self.__seqs[0]
            self.__seqs_new[1] = self.__seqs[1]

    def __mutate_cds(self):
        '''Mutates CDS.'''
        self.__seqs_new[2] = \
            [self.__cod_opt.mutate(prot_seq, dna_seq, 3.0 / len(dna_seq))
             for dna_seq, prot_seq in zip(self.__seqs[2],
                                          self.__prot_seqs.values())]

    def __get_valid_rand_seq(self, length):
        '''Returns a valid random sequence of supplied length.'''
        seq = ''.join([_rand_nuc() for _ in range(0, length)])

        if seq_utils.count_pattern(seq, _INVALID_PATTERN) + \
                seq_utils.count_pattern(seq, '[AGT]TG') == 0:
            return seq

        return self.__get_valid_rand_seq(length)

    def __repr__(self):
        # return '%r' % (self.__dict__)
        cai = [self.__cod_opt.get_cai(prot_seq) for prot_seq in self.__seqs[2]]
        invalid_patterns = [seq_utils.count_pattern(seq, _INVALID_PATTERN)
                            for seq in self.__seqs]
        start_codons = [seq_utils.count_pattern(seq, '[AGT]TG')
                        for seq in self.__seqs]

        return str(cai) + '\t' + str(invalid_patterns) + '\t' + \
            str(start_codons) + '\t' + str(_get_tirs(self.__dgs)) + '\t' + \
            self.__seqs[0] + ' ' + self.__seqs[1] + ' ' + \
            str(len(self.__seqs[2])) + ' '

    def __print__(self):
        return self.__repr__

    def print_sol(self):
        '''Prints the solution.'''
        tirs = _get_tirs(self.__dgs)
        for i, cds in enumerate(self.__seqs[2]):
            print str(self.__dgs[i]) + '\t' + str(tirs[i]) + '\t' + \
                self.__seqs[0] + '\t' + self.__seqs[1] + '\t' + cds + '\t' + \
                self.__seqs[3] + '\t' + \
                ('' if self.__seqs[4] is None else self.__seqs[4])


def _get_tirs(dgs):
    '''Gets the translation initiation rate.'''
    return [_RBS_CALC.calc_expression_level(d_g) for d_g in dgs]


def _replace(sequence, pos, nuc):
    '''Replace nucleotide at pos with nuc.'''
    return sequence[:pos] + nuc + sequence[pos + 1:]


def _rand_nuc():
    '''Returns a random nucleotide.'''
    return random.choice(['A', 'T', 'G', 'C'])


def main(argv):
    '''main method.'''
    print sim_ann.optimise(RBSSolution(argv[5:], argv[1],
                                       len_target=int(argv[2]),
                                       tir_target=float(argv[3])),
                           acceptance=float(argv[4]))


if __name__ == '__main__':
    main(sys.argv)
