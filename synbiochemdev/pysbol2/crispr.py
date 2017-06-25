'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from sbol import ComponentDefinition, Document, Sequence, setHomespace, \
    SO_PROMOTER, SO_CDS


def crispr_test():
    '''Simple example.'''
    setHomespace('http://sys-bio.org')
    doc = Document()

    mkate, _ = _add_comp_def(doc, 'mKate')
    gal4vp16, _ = _add_comp_def(doc, 'gal4vp16')
    pconst, _ = _add_comp_def(doc, 'pconst', SO_PROMOTER, 'atgtaa')
    mkate_cds, _ = _add_comp_def(doc, 'mkate_cds', SO_CDS, 'attcga')
    gal4vp16_cds, _ = _add_comp_def(doc, 'gal4vp16_cds', SO_CDS, 'attcga')
    mkate_gene, mkate_gene_seq = _add_comp_def(doc, 'mkate_gene')
    gal4vp16_gene, gal4vp16_gene_seq = _add_comp_def(doc, 'gal4vp16_gene')

    mkate_gene.assemble([pconst, mkate_cds])
    gal4vp16_gene.assemble([pconst, gal4vp16_cds])

    mkate_gene_seq.assemble()
    print mkate_gene_seq.elements.get()

    gal4vp16_gene_seq.assemble()
    print gal4vp16_gene_seq.elements.get()

    doc.write('crispr.xml')


def _add_comp_def(doc, uri, role=None, seq=None):
    '''Adds ComponentDefinition.'''
    comp_def = ComponentDefinition(uri)
    doc.addComponentDefinition(comp_def)

    if role is not None:
        comp_def.roles.set(role)

    sequence = Sequence(uri, seq) if seq is not None else Sequence(uri)
    doc.addSequence(sequence)

    comp_def.sequence.set(sequence.identity.get())

    return comp_def, sequence


if __name__ == '__main__':
    crispr_test()
