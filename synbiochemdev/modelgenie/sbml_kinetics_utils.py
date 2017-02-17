'''
ModelGenie (c) University of Manchester 2015

ModelGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from StringIO import StringIO
import json

import libsbml

from synbiochemdev.modelgenie import sbml_utils


class KineticsAdder(object):
    ''' Adds (unparameterised) convenience kinetics to SBML Document.'''

    def __init__(self, document, species_to_ignore=None):
        def_species_to_ignore = ['http://identifiers.org/chebi/CHEBI:15377',
                                 'http://identifiers.org/chebi/CHEBI:15378',
                                 'http://identifiers.org/chebi/CHEBI:16526',
                                 'http://identifiers.org/chebi/CHEBI:24636']

        self.__document = document
        self.__species_to_ignore = species_to_ignore \
            if species_to_ignore is not None else def_species_to_ignore
        self.__model = document.getModel()
        self.__units = sbml_utils.get_units()
        self.__formula_to_id = {}

    def add_kinetics(self):
        '''Adds convenience kinetics to all reactions in the document.'''
        for species in self.__model.getListOfSpecies():
            for annotation in sbml_utils.get_annotations(species):
                if annotation[0] in self.__species_to_ignore:
                    species.setConstant(True)
                    species.setBoundaryCondition(True)

        for reaction in self.__model.getListOfReactions():
            self.__add_kinetics(reaction)

    def get_document(self):
        '''Returns the SBML Document.'''
        return self.__document

    def __add_kinetics(self, reaction):
        '''Adds a convenience kinetic law to a given reaction.'''
        is_reversible = reaction.isSetReversible() and reaction.getReversible()

        enzymes = [modifier for modifier in reaction.getListOfModifiers()
                   if modifier.getSboTerm() ==
                   sbml_utils.SBO_TERMS[sbml_utils.ENZYME]]

        if len(enzymes) == 0:
            compartment = self.__model.getListOfCompartments()[0]
            enzyme = self.__add_enzyme(reaction, compartment.getId())
        else:
            enzyme = enzymes[0]

        formula, parameters, num_react, num_prods = \
            self.__get_formula(reaction.getListOfReactants(),
                               reaction.getListOfProducts()
                               if is_reversible else [],
                               enzyme)

        print formula

        func_def = self.__get_func_def(formula,
                                       [parameter[0]
                                        for parameter in parameters],
                                       num_react,
                                       num_prods)

        self.__set_kinetic_law(reaction, func_def, parameters)

    def __add_enzyme(self, reaction, cmpt_id):
        '''Adds enzyme to reaction in given compartment.'''
        sbo_term = sbml_utils.SBO_TERMS[sbml_utils.PROTEIN]
        modifier = sbml_utils.get_species(self.__model, cmpt_id,
                                          sbml_utils.get_unique_id(),
                                          sbo_term)
        ref = reaction.createModifier()
        ref.setSpecies(modifier.getId())
        ref.setSBOTerm(sbml_utils.SBO_TERMS[sbml_utils.ENZYME])
        return modifier

    def __is_constant(self, species_id):
        '''Returns boolean indicating if species is constant.'''
        return self.__model.getSpecies(species_id).getConstant()

    def __get_formula(self, reactants, products, enzyme):
        '''Returns formula, react_terms, prod_terms
        for supplied number of reactants and products.'''
        react_terms = self.__get_terms(reactants, 'S')
        prod_terms = self.__get_terms(products, 'P')
        irreversible = len(prod_terms) == 0

        react_numer_terms = sbml_utils.KCAT_FOR + ' * ' + \
            _get_numer_terms(react_terms)
        react_denom_terms = _get_denom_terms(react_terms)

        prod_numer_terms = '' if irreversible \
            else ' - ( ' + sbml_utils.KCAT_REV + ' * ' + \
            _get_numer_terms(prod_terms) + ' )'
        prod_denom_terms = '' if irreversible \
            else ' + ( ' + _get_denom_terms(prod_terms) + ' ) - 1'

        numer = '( ' + react_numer_terms + prod_numer_terms + ' )'
        denom = '( ' + react_denom_terms + prod_denom_terms + ' )'

        formula = sbml_utils.VOLUME + ' * ' + sbml_utils.ENZYME_CONC + \
            ' * ' + numer + ' / ' + denom

        parameters = []
        parameters.append((sbml_utils.VOLUME,
                           sbml_utils.SBO_TERMS[sbml_utils.VOLUME],
                           enzyme.getCompartment(),
                           None))

        parameters.append((sbml_utils.ENZYME_CONC,
                           sbml_utils.SBO_TERMS[sbml_utils.CONC],
                           enzyme.getId(),
                           None))

        parameters.append((sbml_utils.KCAT_FOR,
                           sbml_utils.SBO_TERMS[sbml_utils.KCAT_FOR],
                           sbml_utils.KCAT_FOR,
                           None))

        if not irreversible:
            parameters.append((sbml_utils.KCAT_REV,
                               sbml_utils.SBO_TERMS[sbml_utils.KCAT_REV],
                               sbml_utils.KCAT_REV,
                               None))

        parameters.extend(_get_parameters(react_terms))
        parameters.extend(_get_parameters(prod_terms))

        return formula, parameters, len(react_terms), len(prod_terms)

    def __get_terms(self, participants, prefix):
        ''''Get list of tuples of
        (id, stoichiometry, parameter_id, sbo_term).'''
        terms = [[participant.getSpecies(),
                  participant.getStoichiometry(),
                  (prefix + str(i + 1),
                   sbml_utils.SBO_TERMS[sbml_utils.CONC]),
                  ('KM_' + prefix + str(i + 1),
                   sbml_utils.SBO_TERMS[sbml_utils.K_M])]
                 for i, participant in enumerate(participants)
                 if not self.__is_constant(participant.getSpecies())]

        return terms

    def __get_func_def(self, formula, parameters, num_reacts, num_prods):
        '''Gets existing or creates new functionDefinition from given
        parameters.'''
        if formula in self.__formula_to_id:
            formula_id = self.__formula_to_id[formula]
            function_definition = \
                self.__model.getFunctionDefinition(formula_id)
        else:
            function_definition = self.__model.createFunctionDefinition()
            function_definition.setId(sbml_utils.get_unique_id())
            function_definition.setName(_get_func_name(num_reacts, num_prods))
            function_definition.setMath(_get_math(formula, parameters))
            self.__formula_to_id[formula] = function_definition.getId()

        return function_definition

    def __set_kinetic_law(self, reaction, function_definition, parameters):
        '''Sets kineticLaw element to reaction.'''

        mathml = '<?xml version="1.0" encoding="UTF-8"?>'
        mathml += '<math xmlns="http://www.w3.org/1998/Math/MathML">'
        mathml += '<apply>'
        mathml += '<ci>' + function_definition.getId() + '</ci>'

        for parameter in parameters:
            mathml += '<ci>' + parameter[2] + '</ci>'

        mathml += '</apply>'
        mathml += '</math>'

        kinetic_law = reaction.createKineticLaw()
        kinetic_law.setMath(libsbml.readMathMLFromString(mathml))

        for parameter in parameters:
            sbo_term = parameter[1]
            if sbo_term != sbml_utils.SBO_TERMS[sbml_utils.CONC] and \
                    sbo_term != sbml_utils.SBO_TERMS[sbml_utils.VOLUME]:
                param = kinetic_law.createParameter()
                param.setId(parameter[2])
                param.setValue(sbml_utils.SBO_TERM_DEFAULT[sbo_term])
                param.setUnits(self.__units[sbo_term])
                param.setSBOTerm(sbo_term)

        return kinetic_law


def parameterise_kinetics(model):
    '''Parameterises all reactions in the document.'''
    kinetic_params = []

    for reaction in model.getListOfReactions():
        kinetic_params.append(_get_kinetic_params(reaction))

    kinetic_params_json = json.dumps(kinetic_params, sort_keys=True, indent=2)
    kinetic_params_json = _get_kinetic_param_values(kinetic_params_json)
    _set_kinetic_param_values(model, kinetic_params_json)


def _get_numer_terms(terms):
    '''Returns numerator terms in the form S1/KM_S1 * S2/KM_S2.'''
    lst = [['( ' + term[2][0] + ' / ' + term[3][0] + ' )'] * int(term[1])
           for term in terms]
    return ' * '.join([item for sublist in lst for item in sublist])


def _get_denom_terms(terms):
    '''Returns denominator terms in the form
    (((S1/KM_S1)^0 + (S1/KM_S1)^1)) * ((S2/KM_S2)^0 + (S2/KM_S2)^1)).'''
    lst = [' + '.join(['( ( ' +
                       term[2][0] +
                       ' / ' +
                       term[3][0] +
                       ' ) ^ ' +
                       str(x) +
                       ' )'
                       for x in range(int(term[1]) + 1)])
           for term in terms]
    return '( ' + ' ) * ( '.join(lst) + ' )'


def _get_parameters(terms):
    '''Gets parameters derived from terms.'''
    lst = [[(term[2][0], term[2][1], term[0], None),
            (term[3][0], term[3][1], 'KM_' + term[0], term[0])]
           for term in terms]
    return [item for sublist in lst for item in sublist]


def _get_func_name(num_reactants, num_products):
    '''Returns a function name based on the number of reactants and
    products.'''
    is_reversible = num_products > 0
    reversible = 'reversible' if is_reversible else 'irreversible'
    name = 'Convenience (' + reversible + '): ' + \
        str(num_reactants) + ' reactants'

    if is_reversible:
        name += ', ' + str(num_products) + ' products'

    return name


def _get_math(formula, parameters):
    '''Returns math element from formula and parameters.'''
    math_elem = '<math xmlns="http://www.w3.org/1998/Math/MathML">'
    param_str = math_elem + '<lambda>'

    for parameter in parameters:
        param_str += '<bvar><ci>' + parameter + '</ci></bvar>'

    mathml = libsbml.writeMathMLToString(libsbml.parseFormula(formula))
    mathml = mathml.replace(math_elem, param_str)
    mathml = mathml.replace('</math>', '</lambda></math>')
    return libsbml.readMathMLFromString(mathml)


def _get_kinetic_params(reaction):
    '''Generates JSON document define kinetic parameters and other details
    required for reaction.'''
    kinetic_params = {}
    parameters = []

    kinetic_params['reaction_id'] = reaction.getId()
    kinetic_params['ec_term'] = \
        sbml_utils.get_typed_annotations(reaction,
                                         'http://identifiers.org/ec-code/')
    kinetic_params['uniprot_id'] = 'unknown'
    kinetic_params['parameters'] = parameters

    kinetic_law = reaction.getKineticLaw()

    for local_parameter in kinetic_law.getListOfParameters():
        parameter = {}

        parameter_id = local_parameter.getId()
        parameter['parameter_id'] = parameter_id

        sbo_term = local_parameter.getSBOTerm()

        if sbo_term == sbml_utils.SBO_TERMS[sbml_utils.K_M]:
            parameter['substrate'] = parameter_id[len('KM_'):]

        parameter['sboTerm'] = sbo_term
        parameter['value'] = local_parameter.getValue()
        parameter['units'] = local_parameter.getUnits()

        parameters.append(parameter)

    return kinetic_params


def _get_kinetic_param_values(kinetic_params_json):
    '''Fills dummy parameters document with real values.'''
    sio = StringIO(kinetic_params_json)
    kinetic_params = json.load(sio)

    for kinetic_param in kinetic_params:
        reaction_id = str(kinetic_param['reaction_id'])

        for parameter in kinetic_param['parameters']:
            parameter['value'] = \
                get_kinetic_param_value(reaction_id,
                                        str(parameter['parameter_id']))

    return json.dumps(kinetic_params, sort_keys=True, indent=2)


def _set_kinetic_param_values(model, kinetic_params_json):
    '''Fills dummy parameters document with real values.'''
    sio = StringIO(kinetic_params_json)
    kinetic_params = json.load(sio)

    for kinetic_param in kinetic_params:
        reaction_id = str(kinetic_param['reaction_id'])
        reaction = model.getReaction(reaction_id)
        kinetic_law = reaction.getKineticLaw()

        for parameter in kinetic_param['parameters']:
            param = kinetic_law.getParameter(str(parameter['parameter_id']))
            param.setValue(parameter['value'])


def get_kinetic_param_value(reaction_id, parameter_id):
    '''Temporary method to add hard coded parameters.'''
    # print '(\'' + reaction_id + '\', \'' + parameter_id + '\'): 1,'

    hard_coded = {
        ('R02013', 'kcat_for'): 3.3e-6,
        ('R02013', 'KM_C00341'): 1,
        ('R02245', 'kcat_for'): 14,
        ('R02245', 'KM_C00002'): 1.1e-3,
        ('R02245', 'KM_C00418'): 3.3e-4,
        ('R00238', 'kcat_for'): 2.1e-3,
        ('R00238', 'kcat_rev'): 8.0e-2,
        ('R00238', 'KM_C00024'): 1,
        ('R00238', 'KM_C00010'): 1,
        ('R00238', 'KM_C00332'): 1,
        ('R01978', 'kcat_for'): 1,
        ('R01978', 'KM_C00024'): 1,
        ('R01978', 'KM_C00332'): 1,
        ('R02082', 'kcat_for'): 1,
        ('R02082', 'KM_C00356'): 1,
        ('R02082', 'KM_C00005'): 1,
        ('R01123', 'kcat_for'): 1,
        ('R01123', 'kcat_rev'): 1,
        ('R01123', 'KM_C00129'): 1,
        ('R01123', 'KM_C00235'): 1,
        ('R01121', 'kcat_for'): 1,
        ('R01121', 'KM_C00002'): 1,
        ('R01121', 'KM_C01143'): 1,
        ('R03245', 'kcat_for'): 1,
        ('R03245', 'KM_C00002'): 1,
        ('R03245', 'KM_C01107'): 1,
        ('R01658', 'kcat_for'): 1,
        ('R01658', 'KM_C00235'): 1,
        ('R01658', 'KM_C00129'): 1
    }

    return hard_coded[(reaction_id, parameter_id)]
