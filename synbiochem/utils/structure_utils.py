'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from matplotlib.colors import LinearSegmentedColormap
import numpy
import pylab
import scipy.spatial
import sys
import tempfile
import urllib

from Bio.PDB.PDBParser import PDBParser


def calc_proximity(pdb_id):
    '''Calculates residue proximities from PDB file.'''
    with tempfile.TemporaryFile() as pdb_file:
        opener = urllib.URLopener()
        opener.retrieve('http://www.rcsb.org/pdb/files/' + pdb_id + '.pdb',
                        pdb_file.name)
        parser = PDBParser()

        chains = [c for c in parser.get_structure(pdb_id,
                                                  pdb_file.name).get_chains()]

        coords = [residue.child_dict['CA'].get_coord()
                  for residue in chains[0] if 'CA' in residue.child_dict]

        return scipy.spatial.distance.cdist(coords, coords, 'euclidean')


def plot(values, plot_filename, plot_format, title, max_value=None):
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


def export_html(values):
    print ','.join(['[' + ','.join([str(v) for v in value.tolist()]) + ']'
                    for value in values])


def main(argv):
    '''main method.'''
    pdb_id = argv[1]
    proximities = calc_proximity(pdb_id)

    # export_html(proximities)

    plot_format = 'png'
    plot_filename = pdb_id + '.' + plot_format
    plot(proximities, plot_filename, plot_format,
         pdb_id + ' proximity plot')


if __name__ == '__main__':
    main(sys.argv)
