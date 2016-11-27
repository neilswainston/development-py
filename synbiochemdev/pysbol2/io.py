'''
synbiochem (c) University of Manchester 2016

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from sbol import Document


def main():
    '''main method.'''
    doc = Document()
    doc.read('sbol.xml')

    for seq in doc.sequences:
        print seq.identity.get()
        print seq.elements.get()

    doc.write('out.xml')


if __name__ == '__main__':
    main()
