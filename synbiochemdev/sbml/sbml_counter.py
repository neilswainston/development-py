'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from collections import defaultdict
import sys

from libsbml import SBMLReader


def main(args):
    '''main method'''
    reader = SBMLReader()
    document = reader.readSBML(args[0])
    model = document.getModel()

    mets = defaultdict(int)
    cmpt_mets = defaultdict(int)
    ids = {}

    for reaction in model.getListOfReactions():
        _count(reaction.getListOfReactants(), mets, cmpt_mets, ids)
        _count(reaction.getListOfProducts(), mets, cmpt_mets, ids)

    for met, count in cmpt_mets.iteritems():
        spec = model.getSpecies(met)
        print '\t'.join([met, spec.getName(), _get_inchi(spec),
                         spec.getCompartment(), str(count)])

    print

    for met, count in mets.iteritems():
        spec = model.getSpecies(ids[met])
        print '\t'.join([met, spec.getName(), _get_inchi(spec), str(count)])


def _count(lst, mets, cmpt_mets, ids):
    '''Counts metabolites.'''
    for ref in lst:
        spec_id = ref.getSpecies()
        cmpt_mets[spec_id] = cmpt_mets.get(spec_id, 0) + 1

        mets[spec_id[:-2]] = mets.get(spec_id[:-2], 0) + 1
        ids[spec_id[:-2]] = spec_id


def _get_inchi(sbase):
    '''Gets annotations.'''
    if sbase.getCVTerms() is not None:
        for cv_term in sbase.getCVTerms():
            for idx in range(cv_term.getNumResources()):
                if 'inchi' in cv_term.getResourceURI(idx):
                    return cv_term.getResourceURI(idx)

    return ''


if __name__ == '__main__':
    main(sys.argv[1:])
