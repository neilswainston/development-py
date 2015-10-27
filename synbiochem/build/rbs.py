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
import synbiochem.optimisation.simulated_annealing as sim_ann
import synbiochem.utils.sequence_utils as seq_utils
import synbiochem.utils.uniprot_utils as uniprot_utils


# Necessary to get constants hidden as class variables in RBS_Calculator:
_RBS_CALC = RBS_Calculator.RBS_Calculator('A', [0, 0])


class RBSSolution(object):
    '''Solution for RBS optimisation.'''

    def __init__(self, protein_ids, taxonomy_id, len_target, tir_target):
        # Check if dg_total or TIR (translation initiation rate) was specified.
        # If TIR, then convert to dg_total.
        self.__dg_target = _RBS_CALC.RT_eff * \
            (_RBS_CALC.logK - math.log(float(tir_target)))
        self.__cod_opt = seq_utils.CodonOptimiser(taxonomy_id)

        self.__prot_seqs = _get_sequences(protein_ids)
        cds = [self.__cod_opt.get_codon_optimised_seq(prot_seq)
               for prot_seq in self.__prot_seqs.values()]

        # Randomly choose an RBS that is a decent starting point,
        # using the first CDS as the upstream sequence:
        (rbs, _) = RBS_MC_Design.GetInitialRBS('', cds[0],
                                               self.__dg_target)

        post_seq_length = 30
        self.__seqs = [_get_valid_rand_seq(max(0, len_target - len(rbs))),
                       rbs,
                       cds,
                       _get_valid_rand_seq(post_seq_length)
                       if self.__prot_seqs > 1 else None]
        self.__dgs = self.__calc_dgs(rbs)
        self.__seqs_new = [None, None, cds, self.__seqs[3]]
        self.__dgs_new = None

    def get_energy(self, dgs=None, cdss=None):
        '''Gets the (simulated annealing) energy.'''
        dgs = self.__dgs if dgs is None else dgs
        cdss = self.__seqs[2] if cdss is None else cdss
        cais = [self.__cod_opt.get_cai(cds) for cds in cdss]
        return sum([abs(d_g - self.__dg_target) for d_g in dgs]) / len(dgs) * \
            (1 - sum(cais) / len(cais)) * \
            sum(_count_invalid_pattern(cdss) +
                _count_pattern(cdss, '[AGT]TG')) / len(cdss)

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
        self.__seqs_new = [None, None, self.__seqs[2], self.__seqs[3]]
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

        if _count_invalid_pattern(pre_seq_new + self.__seqs[1]) + \
                _count_pattern(pre_seq_new + self.__seqs[1], '[AGT]TG') == 0:
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

        if _count_invalid_pattern(pre_seq_new + rbs_new) + \
                _count_pattern(pre_seq_new + rbs_new, '[AGT]TG') == 0:
            self.__seqs_new[0] = pre_seq_new
            self.__seqs_new[1] = rbs_new
        else:
            self.__seqs_new[0] = self.__seqs[0]
            self.__seqs_new[1] = self.__seqs[1]

    def __mutate_cds(self):
        '''Mutates CDS.'''
        self.__seqs_new[2] = \
            [self.__cod_opt.mutate(prot_seq, dna_seq, 1.0 / len(dna_seq))
             for dna_seq, prot_seq in zip(self.__seqs[2],
                                          self.__prot_seqs.values())]

    def __repr__(self):
        # return '%r' % (self.__dict__)
        cai = [self.__cod_opt.get_cai(prot_seq) for prot_seq in self.__seqs[2]]
        invalid_patterns = [_count_invalid_pattern(seq) for seq in self.__seqs]
        start_codons = [_count_pattern(seq, '[AGT]TG') for seq in self.__seqs]

        return str(cai) + '\t' + str(invalid_patterns) + '\t' + \
            str(start_codons) + '\t' + str(_get_tirs(self.__dgs)) + '\t' + \
            self.__seqs[0] + ' ' + self.__seqs[1] + ' ' + \
            str(self.__seqs[2]) + ' ' + self.__seqs[3]

    def __print__(self):
        return self.__repr__


def _get_sequences(protein_ids):
    '''Returns sequences from protein ids, which may be either Uniprot ids,
    or a protein sequence itself.'''
    uniprot_id_pattern = \
        '[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}'

    sequences = {}

    for idx, protein_id in enumerate(protein_ids):
        if re.match(uniprot_id_pattern, protein_id):
            sequences.update(uniprot_utils.get_sequences([protein_id]))
        else:
            sequences[str(idx)] = protein_id

    return sequences


def _get_valid_rand_seq(length):
    '''Returns a valid random sequence of supplied length.'''
    seq = ''.join([_rand_nuc() for _ in range(0, length)])

    if _count_invalid_pattern(seq) + _count_pattern(seq, '[AGT]TG') == 0:
        return seq

    return _get_valid_rand_seq(length)


def _get_tirs(dgs):
    '''Gets the translation initiation rate.'''
    return [_RBS_CALC.calc_expression_level(d_g) for d_g in dgs]


def _replace(sequence, pos, nuc):
    '''Replace nucleotide at pos with nuc.'''
    return sequence[:pos] + nuc + sequence[pos + 1:]


def _rand_nuc():
    '''Returns a random nucleotide.'''
    return random.choice(['A', 'T', 'G', 'C'])


def _count_invalid_pattern(seqs):
    '''Counts invalid patterns in sequence.'''
    max_repeat_nucs = 4
    # Start codons | restriction sites | repeating nucleotides
    pattern = '|'.join(['GGTCTC', 'CACCTGC'] +
                       [x * max_repeat_nucs for x in ['A', 'C', 'G', 'T']])
    return _count_pattern(seqs, pattern)


def _count_pattern(strings, pattern):
    '''Counts pattern in string of list of strings.'''
    if isinstance(strings, str) or isinstance(strings, unicode):
        return len(re.findall(pattern, strings))
    else:
        return [_count_pattern(s, pattern) for s in strings]


def main(argv):
    '''main method.'''
    print sim_ann.optimise(RBSSolution(argv[5:], argv[1],
                                       len_target=int(argv[2]),
                                       tir_target=float(argv[3])),
                           acceptance=float(argv[4]),
                           verbose=False)


if __name__ == '__main__':
    main(sys.argv)
