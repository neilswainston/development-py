'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from itertools import chain
from matplotlib.colors import LinearSegmentedColormap
import pylab
import random
import sys
import tempfile
import urllib

from Bio import SeqUtils
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import numpy
import scipy.spatial
import synbiochemdev.utils.file_utils as file_utils


# KD Hydrophobicity, EIIP, Helix, Sheet, Turn
AMINO_ACID_PROPERTIES = {
    'A': [0.66, 0.336, 0.870, 0.399, 0.326],
    'R': [0.1, 0.707, 0.597, 0.518, 0.569],
    'N': [0.189, 0.123, 0.391, 0.276, 0.9],
    'D': [0.189, 0.9, 0.645, 0.276, 0.855],
    'C': [0.722, 0.625, 0.554, 0.608, 0.719],
    'Q': [0.189, 0.582, 0.691, 0.561, 0.59],
    'E': [0.189, 0.137, 0.9, 0.1, 0.403],
    'G': [0.464, 0.132, 0.213, 0.466, 0.9],
    'H': [0.216, 0.253, 0.691, 0.434, 0.569],
    'I': [0.9, 0.1, 0.589, 0.870, 0.1],
    'L': [0.838, 0.1, 0.814, 0.668, 0.251],
    'K': [0.153, 0.335, 0.755, 0.345, 0.61],
    'M': [0.669, 0.621, 0.828, 0.568, 0.263],
    'F': [0.749, 0.699, 0.683, 0.703, 0.263],
    'P': [0.358, 0.225, 0.1, 0.234, 0.883],
    'S': [0.429, 0.625, 0.323, 0.518, 0.842],
    'T': [0.438, 0.696, 0.391, 0.756, 0.576],
    'W': [0.42, 0.447, 0.622, 0.708, 0.576],
    'Y': [0.384, 0.427, 0.335, 0.746, 0.691],
    'V': [0.873, 0.136, 0.572, 0.9, 0.141]
}


def get_pdb_ids(max_ids=None):
    '''Returns all PDB ids.'''
    url = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx'
    with tempfile.NamedTemporaryFile() as temp:
        urllib.urlretrieve(url, temp.name)
        ids = [line.split()[0] for line in open(temp.name)][2:]

        return ids if max_ids is None \
            else random.sample(ids, min(len(ids), max_ids))


def _get_chain_res_struct(pdb_ids):
    chain_res_struct = {}
    pdb_ids = sorted(pdb_ids)
    pdb_id = pdb_ids.pop(0)
    pdb_id_fasta = '>' + pdb_id
    in_field = False

    with tempfile.NamedTemporaryFile() as temp:
        urllib.urlretrieve('http://www.rcsb.org/pdb/files/ss.txt', temp.name)

        for line in open(temp.name):
            if line.startswith(pdb_id_fasta):
                tokens = line.split('>|:')
                in_field = True
            elif line.startswith('>'):
                if len(pdb_ids):
                    break
                else:
                    pdb_id = pdb_ids.pop(0)
                    pdb_id_fasta = '>' + pdb_id
                    in_field = False
            elif in_field:
                # Do something:
                print line

    return chain_res_struct


def _get_residues(pdb_id):
    '''Gets the residues in a PDB file.'''
    return [SeqUtils.seq1(''.join([residue.get_resname()
                                   for residue in chn
                                   if 'CA' in residue.child_dict]))
            for chn in _get_structure(pdb_id).get_chains()]


def _get_structure(pdb_id):
    '''Returns a PDB structure.'''
    with tempfile.TemporaryFile() as pdb_file:
        opener = urllib.URLopener()
        opener.retrieve('http://www.rcsb.org/pdb/files/' + pdb_id + '.pdb',
                        pdb_file.name)
        parser = PDBParser(QUIET=True)
        return parser.get_structure(pdb_id, pdb_file.name)


def _calc_proximities(pdb_id):
    '''Calculates residue proximities from PDB file.'''
    structure = _get_structure(pdb_id)
    chains = [c for c in structure.get_chains()]

    coords = [[residue.child_dict['CA'].get_coord()
               for residue in chn
               if 'CA' in residue.child_dict]
              for chn in chains]

    return [scipy.spatial.distance.cdist(coord, coord, 'euclidean')
            if len(coord) > 0
            else None
            for coord in coords]


def _get_phi_psi_data(pdb_id):
    '''Gets phi and phi angle data.'''
    builder = PPBuilder()
    return [polypep.get_phi_psi_list()
            for polypep in builder.build_peptides(_get_structure(pdb_id))]


def _plot_proximities(pdb_id):
    '''Plots proximity plot(s).'''
    all_proximities = _calc_proximities(pdb_id)

    plot_format = 'png'

    for idx, proximities in enumerate(all_proximities):
        name = pdb_id + '_' + str(idx + 1)
        _plot(proximities, name + '.' + plot_format, plot_format,
              name + ' proximity plot')


def _get_learning_data(pdb_id):
    '''Gets input/output for deep learning.'''
    all_proximities = _calc_proximities(pdb_id)
    all_phi_psi_data = _get_phi_psi_data(pdb_id)

    # Ensure data consistency in terms of data lengths:
    assert len(all_residues) == len(all_proximities) == \
        len(all_phi_psi_data)
    assert all([len(all_residues[idx]) == len(all_proximities[idx]) ==
                len(all_phi_psi_data[idx])
                for idx in range(len(all_residues))])

    inpt = [list(chain.from_iterable([AMINO_ACID_PROPERTIES[amino_acid]
                                      if amino_acid in AMINO_ACID_PROPERTIES
                                      else [0.0, 0.0, 0.0]
                                      for amino_acid in residues]))
            for residues in all_residues]

    prox_output = [list(chain.from_iterable(proximities))
                   for proximities in all_proximities]

    phi_psi_output = [list(chain.from_iterable(lst))
                      for lst in all_phi_psi_data]

    return pdb_id, all_residues, inpt, prox_output, phi_psi_output


def _sample_learning_data(pdb_id, num_samples, nmer_len):
    '''Samples learning data for deep learning.'''
    pdb_id, all_residues, inpt, prox_output, phi_psi_output = \
        _get_learning_data(pdb_id)

    results = []

    num_aa_props = len(AMINO_ACID_PROPERTIES['A'])

    for _ in range(num_samples):
        chn = int(random.random() * len(all_residues))
        start = int(random.random() * (len(all_residues[chn]) - nmer_len))
        rnge = range(start, start + nmer_len)
        result = [pdb_id,
                  ''.join([all_residues[chn][i] for i in rnge]),
                  inpt[chn][start * num_aa_props:(start + nmer_len) *
                            num_aa_props],
                  [prox_output[chn][i]
                   for i in _get_sub_square_matrix(start, nmer_len,
                                                   len(all_residues[chn]))],
                  phi_psi_output[chn][start * 2:(start + nmer_len) * 2]]

        results.append(result)

    return results


def _plot(values, plot_filename, plot_format, title, max_value=None):
    '''Plots 3d matrix values.'''
    fig = pylab.figure()
    sub_plot = fig.add_subplot(111)

    cmap = LinearSegmentedColormap.from_list(name='name', colors=['g', 'w'],
                                             N=10)

    cax = sub_plot.imshow(values, interpolation='nearest', cmap=cmap)
    cax.set_clim(0.0, max_value)
    sub_plot.set_title(title)

    # Add colorbar, make sure to specify tick locations to match desired tick
    # labels
    min_val = numpy.min(values)
    max_val = numpy.max(values)
    cbar = fig.colorbar(cax, ticks=[min_val, max_val])
    cbar.set_ticks([min_val, max_val])

    pylab.savefig(plot_filename, format=plot_format)


def _get_sub_square_matrix(idx, lngth, size):
    '''Returns square submatrix of length lngth from square matrix of given
    size.'''
    return [((idx + r) * size) + idx + c for r in range(lngth)
            for c in range(lngth)]


def main(argv):
    '''main method.'''
    pdb_ids = get_pdb_ids(int(argv[1]))
    chain_res_struct = _get_chain_res_struct(pdb_ids)

    for pdb_id in pdb_ids:
        try:
            for data in _sample_learning_data(pdb_id, int(argv[2]),
                                              int(argv[3])):
                print '\t'.join([str(datum) for datum in data])
        except AssertionError, err:
            print err


if __name__ == '__main__':
    main(sys.argv)
