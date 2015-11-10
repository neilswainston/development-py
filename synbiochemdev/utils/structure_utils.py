'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from matplotlib.colors import LinearSegmentedColormap
import pylab
import random
import re
import tempfile
import urllib

from Bio import SeqUtils
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import numpy
import scipy.spatial


def get_pdb_ids(max_ids=None):
    '''Returns all PDB ids.'''
    url = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx'
    with tempfile.NamedTemporaryFile() as temp:
        urllib.urlretrieve(url, temp.name)
        ids = [line.split()[0] for line in open(temp.name)][2:]

        return ids if max_ids is None \
            else random.sample(ids, min(len(ids), max_ids))


def get_chain_seq_struct(pdb_ids):
    '''Returns chain, sequence and structure.'''
    chain_seq_struct = {}
    pdb_ids = sorted(pdb_ids)
    pdb_id = pdb_ids.pop(0)
    pdb_id_fasta = '>' + pdb_id
    in_field = False
    tokens = None
    str_data = ''

    with tempfile.NamedTemporaryFile() as temp:
        urllib.urlretrieve('http://www.rcsb.org/pdb/files/ss.txt', temp.name)

        for line in open(temp.name):
            if line.startswith(pdb_id_fasta):
                if in_field:
                    if tokens[:2] not in chain_seq_struct:
                        chain_seq_struct[tokens[:2]] = [None, None]

                    chain_seq_struct[tokens[:2]][0 if tokens[2] == 'sequence'
                                                 else 1] = str_data
                    str_data = ''

                tokens = tuple(re.split('>|:', line.strip())[1:])
                in_field = True

            elif in_field and line.startswith('>'):
                if tokens[:2] not in chain_seq_struct:
                    chain_seq_struct[tokens[:2]] = [None, None]

                chain_seq_struct[tokens[:2]][0 if tokens[2] == 'sequence'
                                             else 1] = str_data
                str_data = ''

                if len(pdb_ids) == 0:
                    break
                else:
                    pdb_id = pdb_ids.pop(0)
                    pdb_id_fasta = '>' + pdb_id
                    in_field = False
                    tokens = None
                    str_data = ''

            elif in_field:
                # Do something:
                str_data += line.strip()

    return chain_seq_struct


def get_sequences(pdb_id):
    '''Gets the sequences in a PDB file.'''
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
        parser = PDBParser(QUIET=True)
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
            if len(coord) > 0
            else None
            for coord in coords]


def get_phi_psi_data(pdb_id):
    '''Gets phi and phi angle data.'''
    builder = PPBuilder()
    return [polypep.get_phi_psi_list()
            for polypep in builder.build_peptides(get_structure(pdb_id))]


def plot_proximities(pdb_id):
    '''Plots proximity plot(s).'''
    all_proximities = calc_proximities(pdb_id)

    plot_format = 'png'

    for idx, proximities in enumerate(all_proximities):
        name = pdb_id + '_' + str(idx + 1)
        _plot(proximities, name + '.' + plot_format, plot_format,
              name + ' proximity plot')


def sample_seq_structs(sample_size, struct_pattern):
    '''Sample sequence and structure data.'''
    sample_seq_structs = []

    while len(sample_seq_structs) < sample_size:
        for chain_seq_struct in get_chain_seq_struct(get_pdb_ids(sample_size)):
            print chain_seq_struct

    return sample_seq_structs[:sample_size]


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
