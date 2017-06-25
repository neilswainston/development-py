'''
ModelGenie (c) University of Manchester 2015

ModelGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import urllib2
import xml.sax

from synbiochemdev.modelgenie import kegg_utils, sbml_utils, \
    sbml_kinetics_utils
import libsbml


def build(pathway_ec_numbers, species_to_ignore=None):
    '''Builds a kinetic model from KEGG pathway entry and EC numbers.'''
    document = kegg_build(pathway_ec_numbers)

    # Annotate
    # Add convenience kinetics
    kinetics_adder = sbml_kinetics_utils.KineticsAdder(document,
                                                       species_to_ignore)
    kinetics_adder.add_kinetics()

    sbml_kinetics_utils.parameterise_kinetics(document.getModel())

    # Check errors: temporary measure until code is mature:
    error_log = sbml_utils.get_error_log(document)
    print error_log.printErrors()

    return document


def kegg_build(pathway_ec_numbers):
    '''Builds a pathway model from KEGG pathway and EC numbers.'''
    document = _init_doc()

    for pathway in pathway_ec_numbers:
        url = 'http://www.kegg.jp/kegg-bin/download?entry=' + pathway + \
            '&format=kgml'

        # Open and parse URL (write bespoke XML parser)
        ec_numbers = pathway_ec_numbers[pathway]
        ec_numbers = ['ec:' + ec_number for ec_number in ec_numbers]
        handler = kegg_utils.KgmlHandler(document, ec_numbers)
        parser = xml.sax.make_parser()
        parser.setContentHandler(handler)
        parser.parse(urllib2.urlopen(url))

    return document


def _init_doc():
    '''Initialises SBML document.'''
    document = libsbml.SBMLDocument(2, 5)
    cmpt_id = 'c'

    model = document.createModel()

    sbml_utils.set_units(model)

    compartment = model.createCompartment()
    compartment.setId(cmpt_id)
    compartment.setSize(1)
    compartment.setUnits('litre')
    sbo_term = sbml_utils.SBO_TERMS[sbml_utils.COMPARTMENT]
    compartment.setSBOTerm(sbo_term)

    return document
