'''
ModelGenie (c) University of Manchester 2015

ModelGenie is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import math
import uuid
import libsbml

KCAT = 'kcat'
KCAT_FOR = 'kcat_for'
KCAT_REV = 'kcat_rev'
K_M = 'KM'
BIOCHEM_REACT = 'BIOCHEM_REACT'
VMAX = 'vmax'
CONC = 'CONC'
COMP_INHIB = 'COMP_INHIB'
NON_COMP_INHIB = 'NON_COMP_INHIB'
SIMPLE_CHEM = 'SIMPLE_CHEM'
PROTEIN = 'PROTEIN'
K_I = 'KI'
KEQ = 'KEQ'
COMPARTMENT = 'COMPARTMENT'
ENZYME = 'ENZYME'
VOLUME = 'volume'
ENZYME_CONC = 'Et'

SBO_TERMS = {KCAT: 25,
             KCAT_FOR: 320,
             KCAT_REV: 321,
             K_M: 27,
             BIOCHEM_REACT: 176,
             VMAX: 186,
             CONC: 196,
             COMP_INHIB: 206,
             NON_COMP_INHIB: 207,
             SIMPLE_CHEM: 247,
             PROTEIN: 252,
             K_I: 261,
             KEQ: 281,
             COMPARTMENT: 290,
             ENZYME: 460,
             VOLUME: 468}

SBO_TERM_DEFAULT = {SBO_TERMS[SIMPLE_CHEM]: 1e-4,
                    SBO_TERMS[PROTEIN]: 1e-6,
                    SBO_TERMS[KCAT]: 10,
                    SBO_TERMS[KCAT_FOR]: 10,
                    SBO_TERMS[KCAT_REV]: 10,
                    SBO_TERMS[K_M]: 1e-4}


def get_species(model, cmpt_id, species_id, sbo_term):
    ''' Gets an existing species, or creates a new one.'''
    if model.getSpecies(species_id) is None:
        species = model.createSpecies()
        species.setId(species_id)
        species.setCompartment(cmpt_id)
        species.setInitialConcentration(SBO_TERM_DEFAULT[sbo_term])
        species.setSBOTerm(sbo_term)

        if sbo_term is SBO_TERMS[PROTEIN]:
            species.setBoundaryCondition(True)
            species.setConstant(True)

        return species

    return None


def add_species_references(model, rctn, cmpt_id, participants):
    ''' Adds SpeciesReferences to a Reaction,
    based on reactant and product definitions.'''
    for participant in participants:
        species_id = participant[0]
        stoich = participant[1]
        sbo_term = SBO_TERMS[SIMPLE_CHEM]
        get_species(model, cmpt_id, species_id, sbo_term)
        ref = rctn.createReactant() if stoich < 0 else rctn.createProduct()
        ref.setSpecies(species_id)
        ref.setStoichiometry(math.fabs(stoich))


def add_annotation(sbase, qual_type, specific_qual_type, resource):
    ''' Create CVTerm and add to SBase.'''
    cv_term = libsbml.CVTerm()
    cv_term.setQualifierType(qual_type)

    if qual_type is libsbml.MODEL_QUALIFIER:
        cv_term.setModelQualifierType(specific_qual_type)
    else:
        cv_term.setBiologicalQualifierType(specific_qual_type)

    cv_term.addResource(resource)
    sbase.setMetaId('_meta_' + sbase.getId())
    sbase.addCVTerm(cv_term)


def get_annotations(sbase):
    '''Gets list of annotations for a given sbase.'''
    annotations = []

    if sbase.isSetAnnotation():
        cv_terms = sbase.getCVTerms()
        for i in range(cv_terms.getSize()):
            cv_term = cv_terms.get(i)
            qual_type = cv_term.getQualifierType()
            specific_qual_type = cv_term.getModelQualifierType() \
                if qual_type == libsbml.MODEL_QUALIFIER \
                else cv_term.getBiologicalQualifierType()

            for j in range(cv_term.getNumResources()):
                resource = cv_term.getResourceURI(j)
                annotations.append((resource, qual_type, specific_qual_type))

    return annotations


def get_typed_annotations(sbase, typ):
    '''Gets list of annotations for a given sbase and type.'''
    return [annotation[0] for annotation in get_annotations(sbase)
            if annotation[0].startswith(typ)]


def set_units(model):
    '''Sets default units on an SBML Model.'''

    # kcat: per second
    kcat_unit_id = 'per_sec'
    unit_def = model.createUnitDefinition()
    unit_def.setId(kcat_unit_id)
    unit_def.setName(kcat_unit_id)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_SECOND)

    # Concentration, K_M, K_I: M
    conc_unit_id = 'M'
    unit_def = model.createUnitDefinition()
    unit_def.setId(conc_unit_id)
    unit_def.setName(conc_unit_id)
    unit = unit_def.createUnit()
    unit.setKind(libsbml.UNIT_KIND_MOLE)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_LITRE)

    # vmax: M per second
    v_unit_id = 'M_per_sec'
    unit_def = model.createUnitDefinition()
    unit_def.setId(v_unit_id)
    unit_def.setName(v_unit_id)
    unit = unit_def.createUnit()
    unit.setKind(libsbml.UNIT_KIND_MOLE)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_LITRE)
    unit = unit_def.createUnit()
    unit.setExponent(-1)
    unit.setKind(libsbml.UNIT_KIND_SECOND)


def get_units():
    '''Sets default units on an SBML Model.'''
    units = {}

    # kcat: per second
    kcat_unit_id = 'per_sec'
    units[SBO_TERMS[KCAT]] = kcat_unit_id
    units[SBO_TERMS[KCAT_FOR]] = kcat_unit_id
    units[SBO_TERMS[KCAT_REV]] = kcat_unit_id

    # Concentration, K_M, K_I: M
    conc_unit_id = 'M'
    units[SBO_TERMS[CONC]] = conc_unit_id
    units[SBO_TERMS[K_M]] = conc_unit_id
    units[SBO_TERMS[K_I]] = conc_unit_id

    # vmax: M per second
    v_unit_id = 'M_per_sec'
    units[SBO_TERMS[VMAX]] = v_unit_id

    # Keq: dimensionless
    units[SBO_TERMS[KEQ]] = 'dimensionless'

    return units


def get_error_log(document):
    '''Checks SBML Document and returns error log.'''
    document.validateSBML()
    document.checkConsistency()
    document.checkInternalConsistency()
    return document.getErrorLog()


def get_unique_id():
    '''Returns unique and SBML-valid id.'''
    return '_' + uuid.uuid4().get_hex().replace('-', '_')
