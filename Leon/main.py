#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
from pathlib import Path
from PIL import Image
from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import shutil
from collections import OrderedDict
from pyrfume import odorants


# We'll work out of the 'patterns_csv' directory, which has glomerular maps + CAS#s bundled together in individual csv files. First flatten the directory:

# In[7]:


saved_path = os.getcwd()
root = '/Users/jcastro/Dropbox/CnG/' # set as needed
patterns_dir = 'Leon/Glomerular archive/patterns_csv' # dir of glomerular maps
path = Path(root + patterns_dir)

# fxn obtained and modified from https://amitd.co/code/python/flatten-a-directory
def flatten(directory):
    for dirpath, _, filenames in os.walk(directory, topdown=False):
        for filename in filenames:
            i = 0
            source = os.path.join(dirpath, filename)
            target = os.path.join(directory, filename)

            while os.path.exists(target):
                i += 1
                file_parts = os.path.splitext(os.path.basename(filename))

                target = os.path.join(
                    directory,
                    file_parts[0] + "_" + str(i) + file_parts[1],
                )

            shutil.move(source, target)

        if dirpath != directory:
            os.rmdir(dirpath)
            
flatten(root+ccc) 


# Walk through the flattened directory, pulling out image data and the corresponding CAS# from the csv files:

# In[8]:


Names = [] # simple file/molecule names (nomenclature idiosyncratic to Leon)
CAS_nums = [] 
Conditions = [] # metadata on experimental conditions (duration, concentration, etc)
Data = [] # list of dfs with the actual data

for file in path.glob('*.csv'):
    df = pd.read_csv(file, index_col=None, header=0)
    # Check to see where data starts
    j=0
    while True:
        if df.iloc[j, 0]=='-100':
            break
        j+=1
    data = df.iloc[j:, :].astype(float)
    data[data<-99] = np.nan
    vals = data.values.ravel()
    assert vals[~np.isnan(vals)].min() > -15
    assert vals[~np.isnan(vals)].max() < 15
    
    name = df.iloc[0, 0]
    cas = df.columns[0]
    condition = df.iloc[1, 0] if j==2 else ''
    
    Names.append(name)
    CAS_nums.append(cas)
    Conditions.append(condition) 
    Data.append(data)


# Clean up the molecule lists, and write csvs of glomerular maps to a new tmp directory. 
# For duplicate maps, just name them as X_replicate_1, X_replicate_2, etc.  

# In[15]:


# first, remove list indices where the CAS# is = 0 (these probably correspond to mixtures, or things like 'coffee', etc)
zeros_idx_list = sorted([ i for i in range(len(CAS_nums)) if CAS_nums[i] == '0' ], reverse=True)

for idx in zeros_idx_list:
    Names.pop(idx)
    CAS_nums.pop(idx)
    Conditions.pop(idx)
    Data.pop(idx)

# cycle through the CAS numbers & identify unique numbers:

lastCAS = CAS_nums[0]
filenames = []

j=0

saved_path = os.getcwd()
root = '/Users/jcastro/Dropbox/CnG/' # set as needed
ccc = 'Leon/Glomerular archive/patterns_csv' # dir of glomerular maps

# make the new directory: 

tmp_path = root + ccc + '_tmp/'
os.mkdir(tmp_path)

for i in range(len(CAS_nums)):
    thisCAS = CAS_nums[i]
    if i>0 and thisCAS == lastCAS: # is the molecule a replicate?
        j = j+1
        name = Names[i] + "_replicate_" + str(j) # ...if so, append a numbered suffix    
    else:
        name = Names[i]
        lastCAS = thisCAS
        j=0
    filenames.append(name)
    Data[i].to_csv(tmp_path + name) # write the glomerular map as a csv (just the image, no metadata)


# For whatever reason, one CAS# (for Methyl 3-aminobenzoate) is just formatted really strangely: 

# In[69]:


# fix a bad CAS#:
idx = CAS_nums.index('10/9/4518')
CAS_nums[idx] = '4518-10-9'


# Fetch CIDs from CAS#s, perform standardization for Pyrfume: 

# In[70]:


molecules = OrderedDict()
molecules = odorants.get_cids(CAS_nums)


# In[71]:


unique_mols = list(molecules.keys())
behavior = []

for cas in unique_mols:
    mol_indices = [i for i in range(len(CAS_nums)) if CAS_nums[i] == str(cas) ] # get all molecules w/ the same CAS#
    behavior.append([filenames[idx] for idx in mol_indices]) # put maps of identical molecules on the same row


# In[82]:


mols = odorants.from_cids(list(molecules.values()))
molecules_df = pd.DataFrame(mols)
molecules_df.set_index('CID', inplace = True)


# Write the 'molecule' and 'behavior' files for pyrfume

# In[124]:


# Molecules
molecules_df.to_csv('molecules.csv')


# In[139]:


import csv

# Behavior
csv_file = open('behavior.csv', 'w')
writer = csv.writer(csv_file)
j = -1

for b in behavior:
    j = j+1
    writer.writerow([list(molecules.values())[j], [b[i] for i in range(len(b))]])





