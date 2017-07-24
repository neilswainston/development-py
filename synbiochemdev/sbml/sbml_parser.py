'''
Created on 24 Jul 2017

@author: neilswainston
'''
from libsbml import SBMLReader


def parse(filename):
    '''Parses SBML file.'''
    values = []
    reader = SBMLReader()
    document = reader.readSBML(filename)
    model = document.getModel()

    for reaction in model.getListOfReactions():
        reac_values = [reaction.getId()]
        cv_terms = reaction.getCVTerms()

        if cv_terms:
            for cv_term in cv_terms:
                for idx in range(cv_term.getNumResources()):
                    uri = cv_term.getResourceURI(idx)

                    if 'ec-code' in uri or 'kegg.reaction' in uri:
                        reac_values.append(uri)

        values.append(reac_values)

    return values


def main():
    '''main method.'''
    for values in parse('recon2_2.xml'):
        print '\t'.join(values)


if __name__ == '__main__':
    main()
