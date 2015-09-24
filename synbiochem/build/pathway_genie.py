'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from synbiochem.build import codon_optimiser, lcr
from synbiochem.utils import uniprot_utils


def build_pathway(target_melt_temp, host_id, uniprot_ids, plasmid_seq=None,
                  shuffle=False, reagent_concs=None):
    '''Returns parts required to build plasmid of pathway using LSR
    technique.'''
    uniprot_seqs = uniprot_utils.get_sequences(uniprot_ids)
    codon_optim = codon_optimiser.CodonOptimiser(host_id)
    parts = _get_parts(uniprot_seqs.values(), codon_optim)
    return lcr.get_bridging_oligos(target_melt_temp,
                                   parts,
                                   plasmid_seq,
                                   shuffle,
                                   reagent_concs)


def _get_parts(protein_seqs, codon_optim):
    '''Gets individual parts.'''
    codon_opt_seqs = codon_optim.optimise(protein_seqs)
    return codon_opt_seqs

print build_pathway(70, '9606', ['P12345', 'P12346'],
                    'tttttttttttttttttttttttttttttttttttttttttttttttttttttttt')
