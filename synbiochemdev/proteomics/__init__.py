'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from collections import OrderedDict, defaultdict
from urllib import urlretrieve
import ast
import copy
import os
import re
import sys

from Bio.SeqUtils.ProtParam import ProteinAnalysis

from pyteomics import mass
from synbiochem.utils import seq_utils
import matplotlib.pyplot as plt


def get_fasta(uniprot_id, filename, variants, digest=None):
    '''Gets a FASTA file representing variants of uniprot id.'''
    uniprot = seq_utils.get_uniprot_values([uniprot_id], ['sequence'])
    get_fasta_from_seq(uniprot[uniprot_id]['Sequence'], filename, variants,
                       digest, uniprot_id)


def get_fasta_from_seq(seq, filename, variants, digest=None, prefix='seq'):
    '''Gets a FASTA file representing variants of sequence.'''
    seqs = []
    apply_variants(seq, [], variants, seqs)

    id_seqs = OrderedDict(
        (prefix + '_' + '_'.join(vals[0]), vals[1]) for vals in seqs)

    id_peps = {seq_id: digest_seq(seq, digest)
               for seq_id, seq in id_seqs.iteritems()}

    unique_id_peps = get_unique_peps(id_peps)
    # print_peps(unique_id_peps)
    draw_spectra(unique_id_peps,  max_charge=3)

    seq_utils.write_fasta(id_seqs, filename)


def apply_variants(seq, var_desc, variants, seqs):
    '''Generates all variants of sequence.'''
    for pos in sorted(variants):
        new_variants = copy.copy(variants)
        new_variants.pop(pos)

        for amino_acid in variants[pos]:
            new_seq = seq[:pos] + amino_acid + seq[pos + 1:]
            new_var_desc = list(var_desc) + \
                [seq[pos] + str(pos + 1) + amino_acid]
            apply_variants(new_seq, new_var_desc, new_variants, seqs)

        return

    seqs.append((var_desc, seq))


def digest_seq(seq, digest):
    '''Digests a sequence according to regex.'''
    return [(pep, ProteinAnalysis(pep).molecular_weight())
            for pep in re.sub(digest, '*',  seq, re.DOTALL).split('*')]


def get_unique_peps(id_peps):
    '''Gets unique peptides.'''
    return {seq_id:
            [pep for pep in peps
             if not all([pep in peps for peps in id_peps.values()])]
            for seq_id, peps in id_peps.iteritems()}


def print_peps(id_peps):
    '''Prints peptides and masses.'''
    for seq_id, peps in id_peps.iteritems():
        print seq_id + '\t' + '\t'.join([str(val)
                                         for pep_mass in peps
                                         for val in pep_mass])


def draw_spectra(id_peps, types=('b', 'y'), max_charge=1):
    '''Draws theoretical fragment spectra.'''
    for seq_id, peps in id_peps.iteritems():
        for pep in peps:
            title = seq_id + ', ' + pep[0] + ' (' + str(pep[1]) + ')'
            plt.figure()
            plt.title(title)
            plt.xlabel('m/z, Th')
            plt.ylabel('Charge')
            plt.ylim(0, max_charge * 1.1)

            frags = fragments(pep[0], types, max_charge)
            plt.bar([frag[0] for frag in frags],
                    [frag[2] for frag in frags],
                    width=0.1,
                    edgecolor=['r' if frag[1] is 'y' else 'b'
                               for frag in frags])
            plt.savefig(title + '.png')


def fragments(peptide, types, max_charge):
    '''The function generates all possible m/z for fragments of types
    and of charges from 1 to maxcharge.'''
    frags = []

    for i in xrange(1, len(peptide) - 1):
        for ion_type in types:
            for charge in xrange(1, max_charge + 1):
                sub_pep = peptide[:i] if ion_type[0] in 'abc' else peptide[i:]
                frags.append((mass.fast_mass(sub_pep, ion_type=ion_type,
                                             charge=charge),
                              ion_type,
                              charge))

    return frags


def main(args):
    '''main method.'''
    variants = ast.literal_eval(args[2])
    get_fasta(
        args[0], args[1], variants, args[3] if len(args) > 3 else None)

if __name__ == '__main__':
    main(sys.argv[1:])
