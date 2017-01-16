'''
Created on 8 Aug 2016

@author: neilswainston
'''
import string

from synbiochem.utils import seq_utils


def _unacceptable(seq, rand_rnge, max_melt_temp=20):
    '''Checks whether seq is acceptable.'''
    id_seqs = {'For': seq, 'Rev': seq_utils.get_rev_comp(seq)}
    results = seq_utils.do_blast(id_seqs, id_seqs, evalue=10, word_size=4)

    for result in results:
        for alignment in result.alignments:
            for hsp in alignment.hsps:
                qu_rnge = range(hsp.query_start, hsp.query_end)
                sb_rnge = range(hsp.sbjct_start, hsp.sbjct_end)

                if hsp.align_length != len(seq) and \
                        len(set(rand_rnge).intersection(qu_rnge)) > 0 and \
                        len(set(rand_rnge).intersection(sb_rnge)) > 0:
                    melt_temp = _get_melting_temp(
                        hsp.query, seq_utils.get_comp(hsp.sbjct))

                    if melt_temp > max_melt_temp:
                        print hsp
                        print melt_temp
                        return True

    return False


def _get_melting_temp(query, subject):
    '''Gets melting temp.'''
    dna1 = ''
    dna2 = ''

    for nucs in zip(query, subject):
        if nucs[0] == '-':
            dna1 = dna1 + seq_utils.get_comp(nucs[1])
        else:
            dna1 = dna1 + nucs[0]

        if nucs[1] == '-':
            dna2 = dna2 + seq_utils.get_comp(nucs[0])
        else:
            dna2 = dna2 + nucs[1]

    return seq_utils.get_melting_temp(dna1, dna2, strict=False)


def main():
    '''main method.'''
    seqs = {}
    seq_id = ''
    seq = ''
    rand_len = 300
    max_repeat_nucl = 3
    invalid_patterns = ['AGATCT', 'CACCTGC', 'CTCGAG', 'GAATTC',
                        'GAGTC([AGCT]{5})', 'GGATCC', 'GGTCTC']

    with open('promoters.txt') as fle:
        for line in fle.read().split('\r'):
            if len(line) > 0 and line[0] == '>':
                if len(seq) > 0:
                    seqs[seq_id] = seq
                seq_id = line[1:].strip()
                seq = ''
            else:
                seq = seq + line.strip().upper()

        seqs[seq_id] = seq

    for seq_id, seq in seqs.iteritems():
        rand_range = range(seq.find('*'), seq.find('*') + rand_len)
        update_seq = string.replace(
            seq, '*', seq_utils.get_random_dna(rand_len, max_repeat_nucl,
                                               invalid_patterns))

        while _unacceptable(update_seq, rand_range):
            update_seq = string.replace(
                seq, '*', seq_utils.get_random_dna(rand_len,
                                                   max_repeat_nucl,
                                                   invalid_patterns))

        seqs[seq_id] = update_seq

    with open('promoters_update.txt', 'w') as fle:
        for seq_id, seq in seqs.iteritems():
            fle.write('>' + seq_id + '\n')
            fle.write(seq + '\n')

if __name__ == '__main__':
    main()
