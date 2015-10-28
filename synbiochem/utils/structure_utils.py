'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from itertools import chain
import math
from matplotlib.colors import LinearSegmentedColormap
import pylab
import random
import scipy.spatial
import sys
import tempfile
import urllib
import urllib2

from Bio import SeqUtils
from Bio.PDB.PDBParser import PDBParser
import numpy


# KD Hydrophobicity, EIIP, Helix, Sheet, Turn
AMINO_ACID_PROPERTIES = {
    'A': [1.8, -0.0667, 32.9, -23.6, -41.6],
    'R': [-4.5, 0.2674, 0, -6.2, -5.1],
    'N': [-3.5, -0.2589, -24.8, -41.6, 44.5],
    'D': [-3.5, 0.4408, 5.8, -41.6, 37.8],
    'C': [2.5, 0.1933, -5.1, 6.8, 17.4],
    'Q': [-3.5, 0.1545, 11.3, 0, -2],
    'E': [-3.5, -0.2463, 36.5, -67.3, -30.1],
    'G': [-0.4, -0.2509, -46.2, -13.9, 44.5],
    'H': [-3.2, -0.1414, 11.3, -18.6, -5.1],
    'I': [4.5, -0.2794, -1, 45.1, -75.5],
    'L': [3.8, -0.2794, 26.2, 15.7, -52.8],
    'K': [-3.9, -0.0679, 19.1, -31.5, 1],
    'M': [1.9, 0.1899, 27.8, 1, -51.1],
    'F': [2.8, 0.26, 10.4, 20.7, -51.1],
    'P': [-1.6, -0.1665, -59.8, -47.8, 41.9],
    'S': [-0.8, 0.1933, -32.9, -6.2, 35.8],
    'T': [-0.7, 0.2572, -24.8, 28.5, -4.1],
    'W': [-0.9, 0.0331, 3, 21.5, -4.1],
    'Y': [-1.3, 0.0148, -31.5, 27, 13.1],
    'V': [4.2, -0.2469, -3, 49.5, -69.3]
}


def get_pdb_ids(scop_id, max_ids=None):
    '''Returns all PDB ids for a given SCOP id.'''
    url = 'http://www.rcsb.org/pdb/results/results.do' + \
        '?startAt=0&resultsperpage=1000000&outformat=text&qrid=' + \
        scop_id

    ids = [line.split()[0] for line in urllib2.urlopen(url)]

    return ids if max_ids is None \
        else random.sample(ids, min(len(ids), max_ids))


def get_residues(pdb_id):
    '''Gets the residues in a PDB file.'''
    return [SeqUtils.seq1(''.join([residue.get_resname()
                                   for residue in chn
                                   if 'CA' in residue.child_dict]))
            for chn in get_structure(pdb_id).get_chains()]


def get_structure(pdb_id):
    '''Returns a PDB structure.'''
    with tempfile.TemporaryFile() as pdb_file:
        opener = urllib.URLopener()
        opener.retrieve('http://www.rcsb.org/pdb/files/' + pdb_id + '.pdb',
                        pdb_file.name)
        parser = PDBParser()
        return parser.get_structure(pdb_id, pdb_file.name)


def calc_proximities(pdb_id):
    '''Calculates residue proximities from PDB file.'''
    structure = get_structure(pdb_id)
    chains = [c for c in structure.get_chains()]

    coords = [[residue.child_dict['CA'].get_coord()
              for residue in chn
              if 'CA' in residue.child_dict]
              for chn in chains]

    return [scipy.spatial.distance.cdist(coord, coord, 'euclidean')
            for coord in coords]


def plot_proximities(pdb_id):
    '''Plots proximity plot(s).'''
    all_proximities = calc_proximities(pdb_id)

    plot_format = 'png'

    for idx, proximities in enumerate(all_proximities):
        name = pdb_id + '_' + str(idx+1)
        _plot(proximities, name + '.' + plot_format, plot_format,
              name + ' proximity plot')


def get_learning_data(pdb_id):
    '''Gets input/output for deep learning.'''
    all_proximities = calc_proximities(pdb_id)
    all_residues = get_residues(pdb_id)

    inpt = [list(chain.from_iterable([AMINO_ACID_PROPERTIES[amino_acid]
                                      if amino_acid in AMINO_ACID_PROPERTIES
                                      else [0.0, 0.0, 0.0]
                                      for amino_acid in residues]))
            for residues in all_residues]

    output = [list(chain.from_iterable(proximities))
              for proximities in all_proximities]

    return pdb_id, all_residues, inpt, output


def sample_learning_data(pdb_id, num_samples, nmer_len):
    pdb_id, all_residues, inpt, output = get_learning_data(pdb_id)

    results = []

    num_aa_props = len(AMINO_ACID_PROPERTIES['A'])

    for _ in range(num_samples):
        chn = int(random.random() * len(all_residues))
        start = int(random.random() * (len(all_residues[chn]) - nmer_len))
        results.append([pdb_id,
                        all_residues[chn][start:start+nmer_len],
                        inpt[chn][start * num_aa_props:(start + nmer_len) *
                                  num_aa_props]])

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


def main(argv):
    '''main method.'''
    for pdb_id in get_pdb_ids(argv[1], argv[2]):
        print '\t'.join([str(x) for x in sample_learning_data(pdb_id,
                                                              int(argv[3]),
                                                              int(argv[4]))])


if __name__ == '__main__':
    main(sys.argv)
