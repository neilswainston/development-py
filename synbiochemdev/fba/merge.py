'''
development-py (c) University of Manchester 2017

development-py is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys
import libsbml


def merge(filenames, out_filename='merge.xml'):
    '''Merge models (from SBML filenames) into a single model.'''
    doc = libsbml.readSBMLFromFile(filenames[0])
    _check(doc)

    for filename in filenames[1:]:
        _merge(doc, libsbml.readSBMLFromFile(filename))

    libsbml.writeSBMLToFile(doc, out_filename)
    _check(doc, True)

    return doc


def _check(doc, recheck=False):
    '''Check SBML document validity.'''
    if recheck:
        doc.checkConsistency()

    doc.printErrors()


def _merge(doc, sub_doc):
    '''Merge submodel into main model.'''
    _check(sub_doc)
    model = doc.getModel()

    sub_model = sub_doc.getModel()

    for sub_species in sub_model.getListOfSpecies():
        _copy_species(model, sub_species)

    for sub_reaction in sub_model.getListOfReactions():
        _copy_reaction(model, sub_reaction)


def _copy_species(model, sub_species):
    '''Copy species.'''
    if not model.getSpecies(sub_species.getId()):
        species = model.createSpecies()
        species.setId(sub_species.getId())
        species.setName(sub_species.getName())
        species.setCompartment(sub_species.getCompartment())
        species.setConstant(False)
        species.setBoundaryCondition(False)


def _copy_reaction(model, sub_reaction):
    '''Copy reaction.'''
    reaction = model.createReaction()
    reaction.setId(sub_reaction.getId())
    reaction.setName(sub_reaction.getName())
    reaction.setReversible(sub_reaction.getReversible())

    lower_bound = _get_bound(model, reaction,
                             _get_lower_bound(sub_reaction))
    upper_bound = _get_bound(model, reaction,
                             _get_upper_bound(sub_reaction),
                             uplow='upper')

    plugin = reaction.getPlugin('fbc')
    plugin.setLowerFluxBound(lower_bound.getId())
    plugin.setUpperFluxBound(upper_bound.getId())

    for sub_reactant in sub_reaction.getListOfReactants():
        reactant = reaction.createReactant()
        reactant.setSpecies(sub_reactant.getSpecies())
        reactant.setStoichiometry(sub_reactant.getStoichiometry())
        reactant.setConstant(True)

    for sub_product in sub_reaction.getListOfProducts():
        product = reaction.createProduct()
        product.setSpecies(sub_product.getSpecies())
        product.setStoichiometry(sub_product.getStoichiometry())
        product.setConstant(True)

    for sub_modifier in sub_reaction.getListOfModifiers():
        modifier = reaction.createModifier()
        modifier.setSpecies(sub_modifier.getSpecies())


def _get_bound(model, reaction, value, units='mmol_per_gDW_per_hr',
               uplow='lower'):
    '''Gets a bound.'''
    bound = model.createParameter()
    bound.setId('_' + uplow + '_bound' + reaction.getId())
    bound.setValue(value)
    bound.setUnits(units)
    return bound


def _get_lower_bound(reaction):
    '''Gets a lower bound.'''
    try:
        return reaction.getLowerFluxBound()
    except AttributeError:
        return -float('Inf') if reaction.getReversible() else 0


def _get_upper_bound(reaction):
    '''Gets an upper bound.'''
    try:
        return reaction.getUpperFluxBound()
    except AttributeError:
        return float('Inf')


def main(args):
    '''main method.'''
    print merge(args)


if __name__ == '__main__':
    main(sys.argv[1:])
