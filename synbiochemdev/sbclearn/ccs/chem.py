'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
import numpy

from rdkit import Chem, DataStructs


def get_fingerprint(smiles, radius=2):
    '''Gets a fingerprint from a SMILES.'''
    mol = Chem.MolFromSmiles(smiles)
    bit_vect = Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius)
    fingerprint = numpy.zeros((1,))
    DataStructs.ConvertToNumpyArray(bit_vect, fingerprint)
    return fingerprint
