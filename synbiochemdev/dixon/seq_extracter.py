'''
synbiochem (c) University of Manchester 2015

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import os
import tempfile
import urllib

from Bio import Seq, SeqIO
from Bio.Alphabet import generic_dna
from numpy import record
from synbiochem.utils import seq_utils


def get_seqs(uniprot_data, outdir, target_len=499):
    '''Gets nucleotide sequences.'''
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    embl_ids = seq_utils.get_uniprot_values(uniprot_data.keys(),
                                            ['database(EMBL)'])

    for uniprot_id, embl_id in embl_ids.iteritems():
        if 'Cross-reference (EMBL)' in embl_id:
            embl_id = embl_id['Cross-reference (EMBL)'][:-1]
            url = 'http://www.ebi.ac.uk/ena/data/view/' + \
                embl_id + '&display=fasta'

            tmpfile = tempfile.NamedTemporaryFile()
            urllib.urlretrieve(url, tmpfile.name)

            with open(tmpfile.name, 'rU') as in_file:
                for record in SeqIO.parse(in_file, 'fasta'):
                    seq = str(record.seq)
                    extn_len = (target_len - len(seq)) / 2
                    forward = uniprot_data[uniprot_id]['strand'] == 'plus'
                    start = max(0, uniprot_data[uniprot_id]['start'
                                                            if forward
                                                            else 'end'] -
                                extn_len)

                    end = min(len(seq),
                              uniprot_data[uniprot_id]['end'
                                                       if forward
                                                       else 'start'] +
                              extn_len)

                    record.id = record.id + '|' + uniprot_id + \
                        '(' + str(start) + ':' + str(end) + ')' + \
                        ('+' if forward else '-')

                    record.description = ''

                    seq = seq[start:end]
                    record.seq = Seq.Seq(seq, generic_dna)

                    with open(os.path.join(outdir, embl_id + '.fasta'),
                              'w') as out_file:
                        SeqIO.write(record, out_file, 'fasta')


def main():
    '''main method.'''
    uniprot_data = {}

    with open('Table S1.txt', 'rU') as in_file:
        for line in in_file:
            tokens = line.split('\t')
            uniprot_data[tokens[1]] = {'start': int(tokens[3]),
                                       'end': int(tokens[4]),
                                       'strand': tokens[5]}

    get_seqs(uniprot_data, 'seqs')


if __name__ == '__main__':
    main()
