
# coding: utf-8

# # Chemical encoding examples

# ## Example 1
# Example of how to encode a chemical compound as input for a neural network.
# * The chemical structure is in SMILES format.
# * Use rdkit to calculate the fingerprint and map into a binary vector.

# In[154]:


import csv

from Bio import SeqIO

from IPython.display import SVG
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
from synbioTools import tensorChem
from synbioTools import tensorSeq


# In[155]:
chems = []
with open('data/solubility/delaney.csv') as f:
    cv = csv.DictReader(f)
    for row in cv:
        chems.append(Chem.MolFromSmiles(row['SMILES']))


# In[156]:


chems[18]


# In[157]:


minPath = 1
maxPath = 5
fpSize = 1024
fp = AllChem.RDKFingerprint(
    chems[18], minPath=1, maxPath=maxPath, fpSize=fpSize)
px = [int(x) for x in list(fp.ToBitString())]
px[0:20]


# We use **synbioTools** to map the chemical into a tensor of shape
# fingerprintSize $\times$ depth. The resulting matrix can be visualized
# using matplotlib.

# In[158]:


depth = 12
fpSize = 20
tc = tensorChem(chems, 20, 4)
tc[0]


# In[159]:


get_ipython().run_line_magic('matplotlib', 'inline')
plt.imshow(tc[0, :, :])
plt.set_cmap('hot')
plt.xlabel('depth')
plt.ylabel('fingerprint')


# ## Example 2
# Example of encoding an amino-acid sequence into a tensor. We use one-hot
# encoding for the amino acid and select desired depth for incuding
# neighboring positions.

# In[160]:


record = list(SeqIO.parse("data/thermostability/l.txt", "fasta"))
seqs = [str(record[i].seq) for i in range(0, len(record))]
seqs[0]


# In[161]:


MAX_SEQ_LENGTH = 10
DEPTH = 5
ts, tss = tensorSeq(seqs, MAX_SEQ_LENGTH, DEPTH)


# In[162]:


get_ipython().run_line_magic('matplotlib', 'inline')
plt.imshow(ts[2])
plt.set_cmap('hot')
plt.xlabel('1-hot kmer')
plt.ylabel('position')


# In[163]:


np.array(aaindex(seqs[0][0:(MAX_SEQ_LENGTH + DEPTH)])


# 1-hot encoding for the first amino acids of sequence 0

# In[164]:


pd.DataFrame(
    np.vstack([np.arange(0, MAX_SEQ_LENGTH + DEPTH),
                np.array(aaindex(seqs[0][0:(MAX_SEQ_LENGTH + DEPTH)]))]
        ), index=['pos', 'seq']
)


# The sequence window is reversed so that any zero-padding appears at the
# beginning of the sequence. In the following example we print the sliding
# window of depth 5 for positions 9,8,7,6:

# In[165]:


ix=np.arange(0, DEPTH * 20, 20)
for a in range(0, 4):
    print(np.where(ts[0, a] == 1) - ix)
