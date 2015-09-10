'''
Created on 10 Sep 2015

@author: neilswainston
'''
import math
import re

NA = 'na'
K = 'k'
TRIS = 'tris'
MG = 'mg'
DNTP = 'dntp'


class MeltingTempCalculator(object):
    '''Class to calculate melting temperatures.'''
    def __init__(self, reagent_concs=None):
        '''See von Ahsen N, Wittwer CT, Schutz E: Oligonucleotide melting
        temperatures under PCR conditions: deoxynucleotide Triphosphate and
        Dimethyl sulfoxide concentrations with comparison to alternative
        empirical formulas. Clin Chem 2001, 47:1956-1961.'''
        na_equivalence_parameter = 3.79

        if reagent_concs is None:
            reagent_concs = {NA: 0.05, K: 0, TRIS: 0, MG: 0.002, DNTP: 1e-7}

        na_equivalence = reagent_concs[NA] + reagent_concs[K] \
            + reagent_concs[TRIS] / 2 \
            + na_equivalence_parameter \
            * math.sqrt(reagent_concs[MG] - reagent_concs[DNTP])

        # See Wetmur J: DNA probes: applications of the principles of nucleic
        # acid hybridization. Crit Rev Biochem Mol Biol 1991, 26:227-259.
        self.base_melting_temp = 81.5 + 16.6 \
            * math.log10(na_equivalence / (1.0 + 0.7 * na_equivalence))

    def get_melting_temp(self, dna1, dna2=None):
        '''Calculates melting temperarure of DNA sequence against its
        complement, or against second DNA sequence.'''
        if dna2 is None:
            dna_gc_content = len(re.findall('G|C', dna1.upper()))
            return self._get_melting_temp(dna_gc_content, 0, len(dna1))
        else:
            # Not implemented...
            return -1

    def _get_melting_temp(self, dna_gc_content, mismatches, length):
        '''Calculates melting temperarure.'''
        gc_increment = 41
        mismatch_decrement = 100
        fixed_decrement = 500
        return self.base_melting_temp + \
            ((gc_increment * dna_gc_content) - fixed_decrement
             - (mismatch_decrement * mismatches)) / length
