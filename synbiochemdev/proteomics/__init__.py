'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from collections import OrderedDict
import ast
import copy
import re
import sys

from synbiochem.utils import seq_utils


def get_fasta(uniprot_id, filename, variants, digest=None):
    '''Gets a FASTA file representing variants of uniprot id.'''
    uniprot = seq_utils.get_uniprot_values([uniprot_id], ['sequence'])
    get_fasta_from_seq(uniprot[uniprot_id]['Sequence'], filename, variants,
                       digest, uniprot_id)


def get_fasta_from_seq(seq, filename, variants, digest=None, prefix='seq'):
    '''Gets a FASTA file representing variants of sequence.'''
    seqs = []
    apply_variants(seq, variants, seqs)

    id_seqs = OrderedDict(
        (prefix + '_' + str(idx), seq) for idx, seq in enumerate(seqs))

    seq_utils.write_fasta(id_seqs, filename)


def apply_variants(seq, variants, seqs):
    '''Generates all variants of sequence.'''
    for pos, amino_acids in variants.iteritems():
        new_variants = copy.copy(variants)
        new_variants.pop(pos)

        for amino_acid in amino_acids:
            new_seq = seq[:pos] + amino_acid + seq[pos + 1:]
            apply_variants(new_seq, new_variants, seqs)

        return

    seqs.append(seq)


def digest_seq(seq, digest):
    '''Digests a sequence according to regex.'''
    return re.sub(digest, '*',  seq, re.DOTALL).split('*')


def main(args):
    '''main method.'''
    variants = ast.literal_eval(args[2])
    get_fasta(
        args[0], args[1], variants, args[3] if len(args) > 3 else None)

if __name__ == '__main__':
    main(sys.argv[1:])
