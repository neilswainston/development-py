
# coding: utf-8

# In[95]:


# import libraries
from __future__ import print_function

import csv
import re

from Bio import Entrez
from Bio import SeqIO
from Bio.PDB.Polypeptide import d1_to_index

import numpy as np
import pandas as pd
from rdkit.Chem import AllChem


# In[38]:
def download_info(patterm, database="Gene"):
    '''download information from relervant database'''

    handle = Entrez.esearch(db=database, term=patterm, retmax=500)
    record = Entrez.read(handle)
    handle.close()

    return record


def download_seq(id_array):
    '''HAVE GeneID to download seqs'''

    result_handle = Entrez.efetch(
        db="nucleotide", rettype="fasta",  id=id_array, retmode="text")
    result = result_handle.read()

    return result


# In[180]:


""" Screen for biosensor sequences """
# List of detectable chemicals from
# https://github.com/brsynth/detectable_metabolites
file_info = 'detectable_metabolites_current.csv'
file_output = 'detectable_metabolites_seq.csv'
Entrez.email = 'dr413677671@gmail.com'

# with open(file_info) as header_input, open(file_output, 'w') as header_output:
#    cw = csv.writer(header_output)
#    cw.writerow( ('Name', 'Inchi', 'Gene', 'Organism', 'NCBI_Gene', 'Name', 'Description', 'Comment', 'Annotation') )
#    cv = csv.DictReader( filter( lambda row: row[0] != '#', header_input ) )
cv = pd.read_csv(file_info, header=0, sep=',')
#    counter=0
record = []

# eliminate unuseful columns
cv.drop(['Reference', 'TypeOfExperiment', 'Comments',
         'TypeOfDetection'], axis=1, inplace=True)
# Change columns names
cv.rename(columns={'DetectedBy': 'Sensor'}, inplace=True)

info = []
id_info = []
seq = []

for index, row in cv.iterrows():
    name = str(row['Name'])
    inchi = row['InChI']
    organism = str(row['Organism'])
    sensor = str(row['Sensor'])
    # Query based on gene name
    # TO DO: curate generic names (riboswitch, TF, etc.)
    if len(sensor) > 0:
        # Query NCBI by gene name
        term = []
        term.append(sensor + '[GN]')
        if organism != 'nan':
            term.append('AND')
            term.append(organism + '[ORGN]')
        term = ' '.join(term)
        try:
            record = download_info(term, database="Gene")
        except:
            continue
        # Fetch the full record or just the summary (should be enough)
        # Things to do:
        # - Select the most convenient format (xml, etc.)
        # - Query by name to avoid false hits
        # - Use Description to double check transcriptional activity
        # - Double-check organims (better through postprocessing)
        if record['IdList']:
            id_array = record['IdList'][0]
            id_info.append(id_array)
            try:
                infotmp, seqtmp = np.squeeze(
                    re.findall(r'>(.*)\n([\s\S]*$)', string))
                info.append(infotmp)
                seq.append(seqtmp)
            except:
                continue
        else:
            info.append("")
            seq.append("")
            id_info.append("")
    else:
        break


# In[91]:


# search 101 sensors only 12 results
info


# In[182]:


# add genes info as new columns to dataframe
cv['info'] = info
cv['id'] = id_info
cv['seq'] = seq

# write csv
file_output = 'detectable_metabolites_seq2.csv'
cv.to_csv(file_output, index=True, sep=',')


# In[181]:


# Count the nonempty seqs num
counter = 0
for i in seq:
    if i:
        counter = counter + 1
print(counter)


# In[ ]:


#             except:
#                 continue

# #                h2 = Entrez.esummary(db='Gene', id=rid, retmode='text')
# #                r2 = [x for x in h2]
#                 r1 = [x.rstrip() for x in h1]
#                 if re.findall(sensor, r1[2]):
#                     try:
#                         h3 = Entrez.efetch(db='sequences', id=rid, rettype='fasta')
#                     except:
#                         continue
#                     out = ( name, inchi, sensor, organism, rid, r1[1], r1[2], r1[5], r1[4] )
#                     print(out)
#                     cw.writerow( out )


# In[ ]:


def readfasta(ffile):
    """ Read fast file, return dictionary """
    record_dict = SeqIO.to_dict(SeqIO.parse(ffile, "fasta"))
    return record_dict


def tensorSeq(seqs, MAX_SEQ_LENGTH, SEQDEPTH, TOKEN_SIZE=20):
    """ Encode an amino acid sequence as a tensor 
    by concatenating one-hot encoding up to desired depth """
    TRAIN_BATCH_SIZE = len(seqs)
    Xs = np.zeros((TRAIN_BATCH_SIZE, MAX_SEQ_LENGTH, SEQDEPTH * TOKEN_SIZE))
    for i in range(0, len(seqs)):
        aaix = aaindex(seqs[i])
        for l in range(0, MAX_SEQ_LENGTH):
            for k in range(0, SEQDEPTH):
                try:
                    Xs[i, l, aaix[l + k] + TOKEN_SIZE * k] = 1
                except:
                    continue
    """ Flip sequences (zero-padding at the start) """
    Xsr = np.flip(Xs, 1)
    return Xsr, Xs


# In[ ]:


def aaindex(seq):
    """ Convert amino acid to numerical index """
    ix = []
    for a in seq:
        if a in d1_to_index:
            ix.append(d1_to_index[a])
    return ix
