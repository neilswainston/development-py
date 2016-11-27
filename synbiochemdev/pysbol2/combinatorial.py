'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from sbol import libsbol, Document, ComponentDefinition, Sequence, \
    SBOL_RESTRICTION_PRECEDES, SBOL_ENCODING_IUPAC


def combinatorial1():
    '''Combinatorial example.'''
    libsbol.setHomespace('http://www.synbiochem.ac.uk')
    doc = Document()

    # SPECIFY POSITIONS #

    # Specify top-level ComponentDefinition and positional Components....
    top = ComponentDefinition('top')
    pos1cd = ComponentDefinition('pos1cd')
    pos2cd = ComponentDefinition('pos2cd')
    pos3cd = ComponentDefinition('pos3cd')
    doc.addComponentDefinition(top)
    doc.addComponentDefinition(pos1cd)
    doc.addComponentDefinition(pos2cd)
    doc.addComponentDefinition(pos3cd)

    pos1c = top.components.create('pos1c')
    pos2c = top.components.create('pos2c')
    pos3c = top.components.create('pos3c')
    pos1c.definition.set(pos1cd.identity.get())
    pos2c.definition.set(pos2cd.identity.get())
    pos3c.definition.set(pos3cd.identity.get())

    # Specify positional SequenceConstraints...
    # pos1 - pos2 - pos3

    pos_sq1 = top.sequenceConstraints.create('pos_sq1')
    pos_sq1.subject.set(pos1c.identity.get())
    pos_sq1.object.set(pos2c.identity.get())
    pos_sq1.restriction.set(SBOL_RESTRICTION_PRECEDES)

    pos_sq2 = top.sequenceConstraints.create('pos_sq2')
    pos_sq2.subject.set(pos2c.identity.get())
    pos_sq2.object.set(pos3c.identity.get())
    pos_sq2.restriction.set(SBOL_RESTRICTION_PRECEDES)

    # SPECIFY 'PARTS' #

    # Specify ComponentDefinitions....
    rbs1cd = ComponentDefinition('rbs1cd')
    rbs2cd = ComponentDefinition('rbs2cd')
    cds1cd = ComponentDefinition('cds1cd')
    cds2cd = ComponentDefinition('cds2cd')

    doc.addComponentDefinition(rbs1cd)
    doc.addComponentDefinition(rbs2cd)
    doc.addComponentDefinition(cds1cd)
    doc.addComponentDefinition(cds2cd)

    # Specify sequences...
    rbs1_seq = Sequence('rbs1cd', 'aaaa', SBOL_ENCODING_IUPAC)
    rbs2_seq = Sequence('rbs2cd', 'cccc', SBOL_ENCODING_IUPAC)
    cds1_seq = Sequence('cds1cd', 'gggg', SBOL_ENCODING_IUPAC)
    cds2_seq = Sequence('cds2cd', 'attcga', SBOL_ENCODING_IUPAC)
    doc.addSequence(rbs1_seq)
    doc.addSequence(rbs2_seq)
    doc.addSequence(cds1_seq)
    doc.addSequence(cds2_seq)

    # Associate sequences with ComponentDefinition...
    rbs1cd.sequence.set(rbs1_seq.identity.get())
    rbs2cd.sequence.set(rbs2_seq.identity.get())
    cds1cd.sequence.set(cds1_seq.identity.get())
    cds2cd.sequence.set(cds2_seq.identity.get())

    # In SBOL Compliant mode, use create methods to add child objects, instead
    # of constructors and add methods (open-world mode)
    rbs1c = top.components.create('rbs1c')
    rbs2c = top.components.create('rbs2c')
    cds1c = top.components.create('cds1c')
    cds2c = top.components.create('cds2c')

#    # Add components...
#    rbs1c = Component('rbs1c')
#    rbs2c = Component('rbs2c')
#    cds1c = Component('cds1c')
#    cds2c = Component('cds2c')
#    top.components.add(rbs1c)
#    top.components.add(rbs2c)
#    top.components.add(cds1c)
#    top.components.add(cds2c)

    # Associate Components to ComponentDefintions...
    rbs1c.definition.set(pos1cd.identity.get())
    rbs2c.definition.set(pos1cd.identity.get())
    cds1c.definition.set(pos2cd.identity.get())
    cds2c.definition.set(pos2cd.identity.get())

    # The following was causing a segmentation fault!!!  Note that rbs1c
    # already belongs to top.components.
    # A child object can't belong to two parent objects.
    # Unfortunately, this is a tricky error to catch...
#    rbs1cd.components.add(rbs1c)
#    rbs2cd.components.add(rbs2c)
#    cds1cd.components.add(cds1c)
#    cds2cd.components.add(cds2c)

    # SPECIFY COMBINATORIAL CONSTRAINTS #

    # Position 1 can be either rbs1c or rbs2...
    combin_pos1_1 = top.sequenceConstraints.create('combin_pos1_1')
    combin_pos1_1.subject.set(pos1c.identity.get())
    combin_pos1_1.object.set(rbs1c.identity.get())
    combin_pos1_1.restriction.set('can_be')

    combin_pos1_2 = top.sequenceConstraints.create('combin_pos1_2')
    combin_pos1_2.subject.set(pos1c.identity.get())
    combin_pos1_2.object.set(rbs2c.identity.get())
    combin_pos1_2.restriction.set('can_be')

#    combin_pos1_1 = SequenceConstraint(
#        'combin_pos1_1', 'pos1c', 'rbs1c', 'can_be')
#    combin_pos1_2 = SequenceConstraint(
#        'combin_pos1_2', 'pos1c', 'rbs2c', 'can_be')
#    top.sequenceConstraints.add(combin_pos1_1)
#    top.sequenceConstraints.add(combin_pos1_2)

    # Position 2 can be either cds1c or cds2c...
    combin_pos2_1 = top.sequenceConstraints.create('combin_pos2_1')
    combin_pos2_1.subject.set(pos2c.identity.get())
    combin_pos2_1.object.set(cds1c.identity.get())
    combin_pos2_1.restriction.set('can_be')

    combin_pos2_2 = top.sequenceConstraints.create('combin_pos2_2')
    combin_pos2_2.subject.set(pos2c.identity.get())
    combin_pos2_2.object.set(cds2c.identity.get())
    combin_pos2_2.restriction.set('can_be')

#    combin_pos2_1 = SequenceConstraint(
#        'combin_pos2_1', 'pos2c', 'cds1c', 'can_be')
#    combin_pos2_2 = SequenceConstraint(
#        'combin_pos2_2', 'pos2c', 'cds2c', 'can_be')
#    top.sequenceConstraints.add(combin_pos2_1)
#    top.sequenceConstraints.add(combin_pos2_2)

    # Position 3 can be either cds1c or cds2c...
    combin_pos3_1 = top.sequenceConstraints.create('combin_pos3_1')
    combin_pos3_1.subject.set(pos3c.identity.get())
    combin_pos3_1.object.set(cds1c.identity.get())
    combin_pos3_1.restriction.set('can_be')

    combin_pos3_2 = top.sequenceConstraints.create('combin_pos3_2')
    combin_pos3_2.subject.set(pos3c.identity.get())
    combin_pos3_2.object.set(cds2c.identity.get())
    combin_pos3_2.restriction.set('can_be')

#    combin_pos3_1 = SequenceConstraint(
#        'combin_pos3_1', 'pos3c', 'cds1c', 'can_be')
#    combin_pos3_2 = SequenceConstraint(
#        'combin_pos3_2', 'pos3c', 'cds2c', 'can_be')
#    top.sequenceConstraints.add(combin_pos3_1)
#    top.sequenceConstraints.add(combin_pos3_2)

    # But the Component at Position 2 cannot be the same as Position 3...
    combin_pos_not_eq = top.sequenceConstraints.create('combin_pos_not_eq')
    combin_pos_not_eq.subject.set(pos2c.identity.get())
    combin_pos_not_eq.object.set(pos3c.identity.get())
    combin_pos_not_eq.restriction.set('not_equal')

#    combin_pos_not_eq = SequenceConstraint(
#        'combin_pos_not_eq', 'pos2c', 'pos2c', 'not_equal')
#    top.sequenceConstraints.add(combin_pos_not_eq)

    for comp_def in doc.componentDefinitions:
        print comp_def.identity.get()

    doc.write('combinatorial.xml')


if __name__ == '__main__':
    combinatorial1()
