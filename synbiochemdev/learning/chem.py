'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=no-member
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
import numpy


def get_fingerprint(smiles, radius=2):
    '''Gets a fingerprint from a SMILES.'''
    mol = Chem.MolFromSmiles(smiles)
    bit_vect = Chem.AllChem.GetMorganFingerprintAsBitVect(mol, radius)
    fingerprint = numpy.zeros((1,))
    DataStructs.ConvertToNumpyArray(bit_vect, fingerprint)
    return fingerprint


def train_rand_for(training_data, values, n_estimators=100, random_state=1123):
    '''Trains a random forest.'''
    # get a random forest classifier:
    rand_for = RandomForestClassifier(n_estimators=n_estimators,
                                      random_state=random_state)
    rand_for.fit(training_data, values)
    return rand_for


def main():
    '''main method.'''
    # generate fingeprints: Morgan fingerprint with radius 2
    fps = [get_fingerprint(smiles) for smiles in ['c1ccccc1',
                                                  'c1ccccc1CC',
                                                  'c1ccncc1',
                                                  'c1ccncc1CC']]

    # train the random forest
    # with the first two molecules being actives (class 1) and
    # the last two being inactives (class 0)
    rand_for = train_rand_for(fps, [1, 1, 0, 0])

    # use the random forest to predict a new molecule
    fingerprint = get_fingerprint('c1ccccc1O')

    print rand_for.predict((fingerprint,))
    print rand_for.predict_proba((fingerprint,))


if __name__ == '__main__':
    main()
