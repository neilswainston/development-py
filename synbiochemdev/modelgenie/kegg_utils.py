'''
ModelGenie (c) University of Manchester 2015

ModelGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import urllib2
import xml.sax

from synbiochemdev.modelgenie import chem_utils, sbml_utils
import libsbml


class KgmlHandler(xml.sax.ContentHandler):
    '''COMMENT'''

    def __init__(self, document, ec_codes):
        self.__ec_codes = ec_codes

        # Instantiate libsbml members:
        self.__document = document
        compartments = document.getModel().getListOfCompartments()
        self.__cmpt_id = compartments[0].getId()
        self.__reaction = None

        self.__reaction_ids = []
        self.__kgml_react_ids = []
        self.__kgml_prod_ids = []

        xml.sax.ContentHandler.__init__(self)

    def get_document(self):
        ''' Returns the SBML Document.'''
        return self.__document

    def startElement(self, name, attrs):
        if name == 'entry' and attrs['name'] in self.__ec_codes:
            self.__reaction_ids.append(attrs['reaction'])
        elif name == 'reaction' and attrs['name'] in self.__reaction_ids:
            # Add reaction:
            level = self.__document.getLevel()
            version = self.__document.getVersion()
            self.__reaction = libsbml.Reaction(level, version)
            reaction_id = attrs['name'].encode('utf8').replace('rn:', '')
            self.__reaction.setId(reaction_id)
            self.__reaction.setReversible(attrs['type'] != 'irreversible')
            sbo_term = sbml_utils.SBO_TERMS[sbml_utils.BIOCHEM_REACT]
            self.__reaction.setSBOTerm(sbo_term)
        elif name == 'substrate' and self.__reaction is not None:
            ref_id = attrs['name'].encode('utf8').replace('cpd:', '')
            self.__kgml_react_ids.append(ref_id)
        elif name == 'product' and self.__reaction is not None:
            ref_id = attrs['name'].encode('utf8').replace('cpd:', '')
            self.__kgml_prod_ids.append(ref_id)

    def endElement(self, name):
        if name == 'reaction' and self.__reaction is not None:
            try:
                # Build reaction:
                model = self.__document.getModel()
                _build_reaction(model, self.__reaction, self.__cmpt_id)

                react_ids = [prt.getSpecies()
                             for prt in self.__reaction.getListOfReactants()]
                prod_ids = [prt.getSpecies()
                            for prt in self.__reaction.getListOfProducts()]

                if all(kgml_reactant_id in prod_ids
                       for kgml_reactant_id in self.__kgml_react_ids) \
                    and all(kgml_product_id in react_ids
                            for kgml_product_id in self.__kgml_prod_ids):
                    # Reaction in KEGG needs to be "flipped":
                    print 'Reaction needs to be flipped'

                elif not all(kgml_reactant_id in react_ids
                             for kgml_reactant_id in self.__kgml_react_ids) \
                    and all(kgml_product_id in prod_ids
                            for kgml_product_id in self.__kgml_prod_ids):
                    # Flag error:
                    msg = 'KEGG equation and KGML definition differ: ' + \
                        self.__reaction.getId()

                    raise ValueError(msg)

                # Add reaction:
                model.addReaction(self.__reaction)

            except ValueError, err:
                print err

            # Clean up:
            self.__kgml_react_ids = []
            self.__kgml_prod_ids = []
            self.__reaction = None

    def endDocument(self):
        xml.sax.ContentHandler.endDocument(self)

        for species in self.__document.getModel().getListOfSpecies():
            species_id = species.getId()
            data = _get_data(species_id)

            species.setName(data['NAME'].split(';')[0])

            resources = []
            resource = 'http://identifiers.org/kegg.compound/' + species_id
            resources.append(resource)
            resources.extend(_get_resources(data['DBLINKS'].split('\n')))

            for resource in resources:
                sbml_utils.add_annotation(species,
                                          libsbml.BIOLOGICAL_QUALIFIER,
                                          libsbml.BQB_IS,
                                          resource)


def _get_data(kegg_id):
    '''Returns KEGG API response in dictionary for a given id.'''
    return __get_data('http://rest.kegg.jp/get/' + kegg_id)


def __get_data(url):
    '''Gets data from the given URL.'''
    data = {}
    key = None

    for line in urllib2.urlopen(url):
        if line[0].isupper():
            key = line.split()[0]
            data[key] = line.replace(key, '').strip()
        else:
            data[key] += ('\n' + line.strip())

    return data


def _build_reaction(model, reaction, cmpt_id):
    '''KGML file doesn't contain the full reaction definition.
    This must be extracted from the KEGG API.'''
    reaction_id = reaction.getId()
    reaction_data = _get_data(reaction_id)

    # Set name:
    reaction.setName(reaction_data['NAME'])

    # Set CV terms:
    resource = 'http://identifiers.org/kegg.reaction/' + reaction_id
    sbml_utils.add_annotation(reaction,
                              libsbml.BIOLOGICAL_QUALIFIER,
                              libsbml.BQB_IS,
                              resource)

    resource = 'http://identifiers.org/ec-code/' + reaction_data['ENZYME']
    sbml_utils.add_annotation(reaction,
                              libsbml.BIOLOGICAL_QUALIFIER,
                              libsbml.BQB_IS_VERSION_OF,
                              resource)

    # Set SpeciesReferences:
    equation = reaction_data['EQUATION']
    participants = chem_utils.parse_equation(equation, '\\<=\\>')
    sbml_utils.add_species_references(model, reaction, cmpt_id, participants)


def _get_resources(dblinks):
    '''Gets resources from given dblinks.'''
    dblinks_map = {'ChEBI': 'http://identifiers.org/chebi/CHEBI:'}
    resources = []

    for dblink in dblinks:
        terms = dblink.split(':')

        if terms[0] in dblinks_map:
            resources.extend([dblinks_map[terms[0]] + link_id
                              for link_id in terms[1].split()])

    return resources
