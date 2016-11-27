'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from sbol import ComponentDefinition, Document, Sequence, setHomespace, \
    SO_PROMOTER, SO_CDS, SO_RBS, SO_TERMINATOR


def example1():
    '''Simple example.'''
    setHomespace("http://sys-bio.org")
    doc = Document()

    gene = ComponentDefinition("BB0001")
    promoter = ComponentDefinition("R0010")
    cds = ComponentDefinition("B0032")
    rbs = ComponentDefinition("E0040")
    terminator = ComponentDefinition("B0012")

    promoter.roles.set(SO_PROMOTER)
    cds.roles.set(SO_CDS)
    rbs.roles.set(SO_RBS)
    terminator.roles.set(SO_TERMINATOR)

    doc.addComponentDefinition(gene)
    doc.addComponentDefinition(promoter)
    doc.addComponentDefinition(cds)
    doc.addComponentDefinition(rbs)
    doc.addComponentDefinition(terminator)

    gene.assemble([promoter, rbs, cds, terminator])

    first = gene.getFirstComponent()
    print first.identity.get()
    last = gene.getLastComponent()
    print last.identity.get()

    # promoter_seq = Sequence("R0010", "ggctgca")
    # rbs_seq = Sequence("B0032", "aattatataaa")
    # cds_seq = Sequence("E0040", "atgtaa")
    # terminator_seq = Sequence("B0012", "attcga")
    # gene_seq = Sequence("BB0001")

    # doc.addSequence(promoter_seq)
    # doc.addSequence(cds_seq)
    # doc.addSequence(rbs_seq)
    # doc.addSequence(terminator_seq)
    # doc.addSequence(gene_seq)

    # promoter.sequence.set(promoter_seq.identity.get())
    # cds.sequence.set(cds_seq.identity.get())
    # rbs.sequence.set(rbs_seq.identity.get())
    # terminator.sequence.set(terminator_seq.identity.get())
    # gene.sequence.set(gene_seq.identity.get())

    # gene_seq.assemble()
    # print gene_seq.elements.get()

    doc.write("gene_cassette.xml")

    # Round trip...
    doc.write('gene_cassette.xml')

    new_doc = Document()
    new_doc.read('gene_cassette.xml')
    new_doc.write('gene_cassette_out.xml')


if __name__ == '__main__':
    example1()
