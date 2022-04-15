#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from pathlib import Path
from IPython.display import display
import numpy as np
import os
from collections import OrderedDict
from pyrfume import odorants


# ## Grab the data from the glomerular archive, and put into a DF

# In[2]:


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


# In[3]:


Names = [] # simple file/molecule names (nomenclature idiosyncratic to Leon)
CAS_nums = [] # CAS numbers (given for each file)
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


# Fix a couple pathological cases, and deal separately w/ 1) monomolecular odorants w/ a CAS vs. 2) complex odorants w/ no CAS (i.e. 'coffee', 'apple')

# In[4]:


# fix a bad CAS#:
idx = CAS_nums.index('10/9/4518')
CAS_nums[idx] = '4518-10-9'


# In[5]:


d = {'Stimulus' : Names, 'CAS numbers': CAS_nums, 'Conditions': Conditions, 'Data': Data}
df = pd.DataFrame(data=d)
df = df.sort_values(by=['CAS numbers', 'Stimulus'])
CAS_numbers = df['CAS numbers'].to_numpy()

# pull out the problem CIDs (cases like 'strawberry', where there will be no CID)
noCID = df[CAS_numbers == '0']
df = df[CAS_numbers != '0']


# Create unique file names for all glomerular images. In cases where there is only one presentation of a given stimulus, we call that file odorantX_0. In cases where there are multiple presentations, we start numbering w/ '1', for example: odorantY_1, odorantY_2, odorantY_3, etc... 

# In[6]:


# create a unique file string for each molecule, grouped by CAS:
grouped = df.groupby((df['CAS numbers'].shift() != df['CAS numbers']).cumsum())

replicates = []

for _,b in grouped:
    x = b['Stimulus'].to_numpy()
    if len(x) == 1:
        replicates.append('csvs/'+ str(x[0]) + '_0.csv')
    else:
        for i in range(len(x)):
            replicates.append('csvs/' + (str(x[i]) + '_' + str(i+1)) + '.csv')

# make the filenames a column in the df, so they follow the data/conditions:            
df['File'] = replicates


# In[7]:


root = '/Users/jcastro/Dropbox/CnG/' # set as needed
local = 'Leon/Glomerular archive/NewTest/'
path = root + local
print(type(path))


# Generating filenames as above, but for the data where there is no CID:

# In[8]:


# same naming convention as w/ the CID-indexed data above
filepaths = ['csvs/' + s + '_0.csv' for s in noCID['Stimulus']]
noCID['File'] = filepaths

# there are only a small number of cases of duplicates, so we'll just hand-fix these
noCID['File'][13] = 'csvs/banana_1.csv'
noCID['File'][79] = 'csvs/banana_2.csv'
noCID['File'][293] = 'csvs/banana_3.csv'

noCID['File'][201] = 'csvs/fujiapple_1.csv'
noCID['File'][213] = 'csvs/fujiapple_2.csv'


# In[9]:


noCID.head(10)


# ## Write csv files (glomerular maps) to disk:

# In[10]:


# write the CAS-indexed data to csv:
for ind in df.index:
    thisdf = df['Data'][ind]
    filepath = root + local + df['File'][ind]
    thisdf.to_csv(filepath,index=False,header=False)

# write the orphans (no CAS#: 'banana', 'wheatbran', etc) to csv:
for ind in noCID.index:
    thisdf = noCID['Data'][ind]
    filepath = root + local + noCID['File'][ind]
    thisdf.to_csv(filepath,index=False,header=False)


# ## Generate standardized dictionary of molecules
# 

# In[11]:


molecules = OrderedDict()
molecules = odorants.get_cids(df['CAS numbers'])


# ## Standardize per pyrfume specifications: 
# 
# Identifier and imaging datasets

# In[12]:


# Make a column of CIDs for the dataframe:
CIDs = [molecules[x] for x in df['CAS numbers']]
df['CIDs'] = CIDs

# dummy cid values for the non-monomolecular odorants
dummy_cids = [-1] * 3 + [-1 * x for x in range(2,7)] + [-7] * 2 + [-1 * x for x in range(8,25)]
noCID['CIDs'] = dummy_cids

# generate the 'identifiers' and 'imaging' datasets
mergedData = pd.concat([noCID, df]).sort_values(by='CIDs')
identifiers = mergedData[['Stimulus', 'CIDs', 'Conditions']]
imaging = mergedData[['Stimulus', 'File']]

# use the names from Leon as the indices
identifiers.set_index('Stimulus', inplace=True)
imaging.set_index('Stimulus', inplace=True)


# Molecules dataset

# In[13]:


# standardized molecules library
# Case 1: molecules w/ CIDs:

mols = odorants.from_cids(list(molecules.values()))
molecules_df = pd.DataFrame(mols)
molecules_df.set_index('CID', inplace = True)

# Case 2: molecules w/o CIDs:
CID_orphans = pd.DataFrame({'CID': dummy_cids,'MolecularWeight':[np.nan]*len(dummy_cids), 'IsomericSMILES':[np.nan]*len(dummy_cids), 'IUPACName':[np.nan]*len(dummy_cids), 'name':noCID['Stimulus']})
CID_orphans.set_index('CID', inplace=True)

# Merged molecule library
all_molecules = pd.concat([molecules_df, CID_orphans]).sort_values(by='CID')


# ## Write pyrfume files to disk:

# In[14]:


all_molecules.to_csv('molecules.csv')
imaging.to_csv('imaging.csv')
identifiers.to_csv('identifiers.csv')

