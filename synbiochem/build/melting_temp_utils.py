'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math
import re

NA = 'NA'
K = 'K'
TRIS = 'TRIS'
MG = 'MG'
DNTP = 'DNTP'


def get_melting_temp(dna1, dna2=None, reagent_concs=None):
    '''Calculates melting temperarure of DNA sequence against its
    complement, or against second DNA sequence.

    See von Ahsen N, Wittwer CT, Schutz E: Oligonucleotide melting
    temperatures under PCR conditions: deoxynucleotide Triphosphate and
    Dimethyl sulfoxide concentrations with comparison to alternative
    empirical formulas. Clin Chem 2001, 47:1956-1961.'''

    base_melt_temp = _get_base_melting_temp(reagent_concs)

    if dna2 is None:
        dna_gc_content = len(re.findall('G|C', dna1.upper()))
        return _get_melting_temp(dna_gc_content, 0, len(dna1), base_melt_temp)
    else:
        # Not implemented...
        return -1


def _get_base_melting_temp(reagent_concs=None):
    '''Gets the base melting temperature based on reagent concentrations.'''
    na_equivalence_parameter = 3.79

    if reagent_concs is None:
        reagent_concs = {NA: 0.05, K: 0, TRIS: 0, MG: 0.002, DNTP: 1e-7}

    na_equiv = reagent_concs[NA] + reagent_concs[K] \
        + reagent_concs[TRIS] / 2 \
        + na_equivalence_parameter \
        * math.sqrt(reagent_concs[MG] - reagent_concs[DNTP])

    # See Wetmur J: DNA probes: applications of the principles of nucleic
    # acid hybridization. Crit Rev Biochem Mol Biol 1991, 26:227-259.
    return 81.5 + 16.6 * math.log10(na_equiv / (1.0 + 0.7 * na_equiv))


def _get_melting_temp(dna_gc_content, mismatches, length, base_melting_temp):
    '''Calculates melting temperature.'''
    gc_increment = 41
    mismatch_decrement = 100
    fixed_decrement = 500
    return base_melting_temp + \
        ((gc_increment * dna_gc_content) - fixed_decrement -
         (mismatch_decrement * mismatches)) / float(length)
