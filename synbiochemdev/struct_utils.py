'''
GeneGenie2 (c) GeneGenie Bioinformatics Limited 2019

All rights reserved.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=wrong-import-order
from itertools import product
import os
import random
import re
from urllib.request import urlopen

from Bio import SeqUtils
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

import numpy as np
import scipy.spatial.distance as dist


_DIR = 'structure_utils'
_PDB_DIR = 'pdb'
_SEC_STRUCT = 'EBHIGTSL '


def get_pdb_ids(max_ids=None, local_only=False):
    '''Returns PDB ids.'''
    if local_only:
        # Returns local PDB ids.
        pdb_dir = os.path.join(os.path.expanduser('~'), _DIR, _PDB_DIR)
        ids = [filename[:-4].upper()
               for _, _, files in os.walk(pdb_dir)
               for filename in files if filename.endswith('.pdb')]
    else:
        # Returns all PDB ids.
        source_url = 'http://www.uniprot.org/uniprot/?query=database:pdb' + \
            '&format=tab&columns=id,database(PDB)'

        with urlopen(source_url) as fle:
            ids = [x for line in fle
                   for x in line.split()[1].split(';')
                   if x and x != 'Cross-reference']

    return ids if max_ids is None \
        else random.sample(ids, min(len(ids), max_ids))


def get_seq_structs(pdb_ids=None):
    '''Returns sequence and structure.'''
    seq_structs = {}
    pdb_ids = sorted(pdb_ids) if pdb_ids is not None else None
    in_field = False
    tokens = ()
    str_data = ''

    source_url = 'http://www.rcsb.org/pdb/files/ss.txt'

    with urlopen(source_url) as fle:
        for line in fle:
            if line.startswith('>'):
                pdb_id = re.search('(?<=\\>)[^:]*', line).group(0)

                if pdb_ids is None or pdb_id in pdb_ids:
                    if in_field:
                        if tokens[:2] not in seq_structs:
                            seq_structs[tokens[:2]] = [None, None]

                        seq_structs[tokens[:2]][0 if tokens[2] == 'sequence'
                                                else 1] = str_data
                        str_data = ''

                    tokens = tuple(re.split('>|:', line.strip())[1:])
                    in_field = True

                elif in_field:
                    if tokens[:2] not in seq_structs:
                        seq_structs[tokens[:2]] = [None, None]

                    seq_structs[tokens[:2]][0 if tokens[2] == 'sequence'
                                            else 1] = str_data
                    str_data = ''

                    in_field = False
                    tokens = ()
                    str_data = ''

            elif in_field:
                str_data += line[:-1]

    return {key: value for key, value in seq_structs.items()
            if all(val is not None for val in value)}


def get_pep_structs(struct_set, length):
    '''Returns sequence and structure.'''
    for struct in struct_set:
        assert struct in _SEC_STRUCT

    return [pep for struct in struct_set for pep in _get_peps(struct, length)]


def get_sequences(pdb_id, chain=None):
    '''Gets the sequences in a PDB file.'''
    return [SeqUtils.seq1(''.join([residue.get_resname()
                                   for residue in chn
                                   if 'CA' in residue.child_dict]))
            for chn in get_structure(pdb_id).get_chains()
            if chain is None or chain == chn.get_id()]


def get_structure(pdb_id):
    '''Returns a PDB structure.'''
    source_url = 'http://www.rcsb.org/pdb/files/' + pdb_id + '.pdb'

    with urlopen(source_url) as fle:
        parser = PDBParser(QUIET=True)
        return parser.get_structure(pdb_id, fle.name)


def get_structure_from_file(filename):
    '''Returns a PDB structure from file.'''
    pdb_id, _ = os.path.splitext(os.path.basename(filename))

    with open(filename) as fle:
        parser = PDBParser(QUIET=True)
        return pdb_id, parser.get_structure(pdb_id, fle.name)


def calc_proximities(structure):
    '''Calculates residue proximities from PDB file.'''
    chain_prox = {}

    for chn in structure.get_chains():
        residues = [SeqUtils.seq1(residue.resname)
                    for residue in chn
                    if 'CA' in residue.child_dict]

        if residues:
            coords = np.array(
                [np.array([atom.get_coord()
                           for atom in residue.child_dict.values()])
                 for residue in chn
                 if 'CA' in residue.child_dict])

            dists = np.array(
                [np.min(dist.cdist(res1, res2))
                 for res1, res2 in product(coords, repeat=2)]
            ).reshape(len(coords), len(coords))

            chain_prox[chn.id] = {'ids': residues,
                                  'seq': ''.join(residues),
                                  'dists': dists}

    return chain_prox


def get_phi_psi_data(pdb_id, chain=None):
    '''Gets phi and phi angle data.'''
    builder = PPBuilder()
    return [polypep.get_phi_psi_list()
            for model in get_structure(pdb_id)
            for chn in model
            if chain is None or chain == chn.get_id()
            for polypep in builder.build_peptides(chn)]


def _get_peps(struct, length):
    '''Gets n-mers of length, matching a given structure.'''
    directory = os.path.join(os.path.expanduser('~'), _DIR, str(length))
    filename = _get_struct_filename(directory, struct)

    if not os.path.isfile(filename):
        _write_peps(length, directory)

    with open(filename) as fle:
        return [line.split()[0] for line in fle]


def _write_peps(length, directory):
    '''Writes n-mers of length, matching a given structure, to file.'''
    assert length % 2 == 1

    pep_structs = {}

    for seq_struct in get_seq_structs().values():
        seq = seq_struct[0]
        struct = seq_struct[1]

        for i in range(len(seq) - length + 1):
            pep_seq = seq[i:i + length]
            pep_struct = struct[i:i + length]

            if pep_seq not in pep_structs:
                pep_structs[pep_seq] = set([pep_struct])
            else:
                pep_structs[pep_seq].add(pep_struct)

    classified_peps = _get_classified_peps(pep_structs, length / 2)
    _write_classified_peps(classified_peps, directory)


def _get_classified_peps(pep_structs, middle_index):
    '''Gets peptides classified by middle residue secondary structure,
    e.g. H, S, B, etc.'''
    classified_peps = {key: [] for key in list(_SEC_STRUCT) + ['unknown']}

    for pep, structs in pep_structs.items():
        structs = list(structs)
        clss = 'unknown' if len(structs) > 1 else structs[0][middle_index]
        classified_peps[clss].append([pep, structs])

    return classified_peps


def _write_classified_peps(classified_peps, directory):
    '''Writes classified peptides to file.'''
    if not os.path.exists(directory):
        os.makedirs(directory)

    for clss, peps in classified_peps.items():
        with open(_get_struct_filename(directory, clss), 'w+') as outfile:
            for pep in peps:
                outfile.write('\t'.join([str(val) for val in pep]) + '\n')


def _get_struct_filename(directory, struct):
    '''Gets filename for a given secondary structure class.'''
    # struct = struct if struct is not ' ' else 'blank'
    return os.path.join(directory, struct + '.xls')
