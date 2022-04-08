#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import scipy.io as sio
import pandas as pd
from pyrfume import odorants
import numpy as np
from collections import OrderedDict
import matplotlib.pyplot as plt


# # I. Glomerular data 
# 
# ### Grabbing the glomerular data

# In[3]:


saved_path = os.getcwd()
root = '/Users/jcastro/Documents/GitHub/pyrfume-data/chae_2019/' # set as needed
glom = 'mosaic_representations/glomerular_data/glomerular_responses/glomerular_responses__main'

if not os.path.exists(root + 'pyrfumeData'):
    os.mkdir(root + 'pyrfumeData') # new directory to store csvs for Pyrfume

os.chdir(root + glom)

glom_responses_left = OrderedDict()
glom_responses_right = OrderedDict()
ROIs_left = OrderedDict()
ROIs_right = OrderedDict()


# walk through the directories and grab the glomerular response data
for root, dirs, files in os.walk(root + glom):
    for filename in files:
        res = root.split('/')
        hemibulb, animal = res[-1], res[-2]
        
        if filename.endswith('responses.mat'):
            
            thedata = sio.loadmat(root + '/'+ filename) 
            
            if hemibulb == 'left_hemibulb':
                glom_responses_left[animal] = thedata['left_hemibulb_ROI_responses']
                
            if hemibulb == 'right_hemibulb':
                glom_responses_right[animal] = thedata['right_hemibulb_ROI_responses']
    
        if filename.endswith('descriptors.mat'):
            
            thedata = sio.loadmat(root + '/' + filename)
            
            if hemibulb == 'left_hemibulb':
                ROIs_left[animal] = thedata['left_hemibulb_ROI_descriptors']
                
            if hemibulb == 'right_hemibulb':
                ROIs_right[animal] = thedata['right_hemibulb_ROI_descriptors']

# grab the molecule list, which is up a directory

root = '/Users/jcastro/Documents/GitHub/pyrfume-data/chae_2019/' 
path = root + 'mosaic_representations/supplementary table 1 odors with cas and cid.xlsx'

glom_CIDs = pd.read_excel(path, header = 0, sheet_name='Glomerular odor panel',
                              usecols = 'C', squeeze=True)

# putting  a dummy value in for the nan, so that the function odorants.from_cids doesn't throw an error
glom_CIDs = glom_CIDs.replace(np.nan, 123456).values.tolist()

#values.tolist()

glom_CIDs = list(map(int,glom_CIDs))
os.chdir(saved_path)


# Generate the CIDs for the molecules:
# ***

# In[4]:


molecules = OrderedDict()
mols = odorants.from_cids(glom_CIDs)

# address the pathological row 'coffee' (row #40), which has no PCID. Obviously. 
mols[40]['CID'] = ''
mols[40]['MolecularWeight'] =''
mols[40]['IsomericSMILES'] =''
mols[40]['IUPACName'] = ''
mols[40]['name'] = 'coffee'

mols_df = pd.DataFrame(mols)
mols_df.set_index('CID', inplace = True)


# In[5]:


mols_df.head()


# ### Write the glomerular data to disk

# In[9]:


# write the molecules to csv:

root = '/Users/jcastro/Documents/GitHub/pyrfume-data/chae_2019/PyrfumeData/' # set as needed
mols_df.to_csv(root + 'molecules_glom_57.csv')

# write the glomerular data to csv

for n in range(1,6):
    animalNum = 'animal_#' + str(n)
    data_left = pd.DataFrame(data = np.array(glom_responses_left[animalNum]))
    data_right = pd.DataFrame(data = np.array(glom_responses_right[animalNum]))
    
    filename = animalNum + '_glom.csv'
    data_left.to_csv(header=False, index=False, path_or_buf = root + 'left_' + filename)
    data_right.to_csv(header=False, index=False, path_or_buf = root + 'right_' + filename)


# # II Mitral/Tufted data
# 
# ## Grabbing the mitral/tufted data

# In[10]:


saved_path = os.getcwd()
root = '/Users/jcastro/Documents/GitHub/pyrfume-data/chae_2019/' 
mitral = 'mosaic_representations/MT_cell/'

#mitral_cells/'
os.chdir(root + mitral)

mitral_33 = OrderedDict()
mitral_55 = OrderedDict()
tufted = OrderedDict()

# walk through the directories and grab the glomerular response data
for root, dirs, files in os.walk(root + mitral):
    for filename in files:
        res = root.split('/')
        experiment_type = res[-1] # tufted, or mitral1 (33 odorants), or mitral2 (55 odorants)
        
        if filename.endswith('.mat'):
            namestring = filename.split('.')[-2]
            thedata = sio.loadmat(root + '/'+ filename) 
            
            if experiment_type == 'tufted_cells':
                tufted[namestring] = thedata['R']
            if experiment_type == '33_odors':
                mitral_33[namestring] = thedata['R']
            if experiment_type == '55_odors':
                mitral_55[namestring] = thedata['R']

# beans = pd.read_excel('MT_cells.xlsx', sheet_name = 'odor list_33 odors', header=0, usecols = 'B').to_numpy()
mitral33_mol_list = pd.read_excel('MT_cells.xlsx', sheet_name='odor list_33 odors', header=0, usecols = 'B').to_numpy()
mitral55_mol_list = pd.read_excel('MT_cells.xlsx', sheet_name='odor list_55 odors', header=0, usecols = 'B').to_numpy()

# flatten:
mitral_33_mol_list = [item[0] for item in mitral33_mol_list]
mitral_55_mol_list = [item[0] for item in mitral55_mol_list]

os.chdir(saved_path)


# In[11]:


# get cids from molecule names, for the panel of 33 odorants: 
mitral_33_cids = odorants.get_cids(mitral_33_mol_list)


# In[12]:


# get cids from molecule names, for the panel of 55 odorants:
mitral_55_cids = odorants.get_cids(mitral_55_mol_list)


# ### Addressing a couple pathological cases where CIDs weren't found: 
# 

# In[13]:


# hard-coding CIDS for pathological cases:

mitral_33_cids['ethyl 3-mercapto propionate'] = 21625
mitral_33_cids['4-siopropyl benzaldehyde'] = 325

mitral_55_cids['3-mercapto propionate'] = 21625
mitral_55_cids['fenchone (-)'] = 82229
mitral_55_cids['citral cis+trans'] = 638011
mitral_55_cids['3-acetal 2,5 dimethyl furan'] = 61527
mitral_55_cids['phenyl ethyl acetate'] = 7654
mitral_55_cids['carvyl acetate (-)'] = 735
mitral_55_cids['2-hexenol'] = 5318042


# ### Write the molecule files for the mitral/tufted data:

# In[16]:


# DataFrames for the various odorant panels. 
# Note that the tufted cell panel is the same as the mitral 55 panel. 

mols_33_mitral = odorants.from_cids(list(mitral_33_cids.values()))
mols_33_mitral_df = pd.DataFrame(mols_33_mitral)
mols_33_mitral_df.set_index('CID', inplace = True)

mols_55_mitral = odorants.from_cids(list(mitral_55_cids.values()))
mols_55_mitral_df = pd.DataFrame(mols_55_mitral)
mols_55_mitral_df.set_index('CID', inplace = True)

mols_tufted_df = mols_55_mitral_df # making this just to avoid ambiguity


# In[22]:


# Write the data: 
root = '/Users/jcastro/Documents/GitHub/pyrfume-data/chae_2019/pyrfumeData/' 

filename = 'molecules_mitral_55.csv'
mols_55_mitral_df.to_csv(root + filename)
filename = 'molecules_mitral_33.csv'
mols_33_mitral_df.to_csv(root + filename)
filename = 'molecules_tufted_55.csv'
mols_tufted_df.to_csv(root + filename)


# In[15]:


root = '/Users/jcastro/Documents/GitHub/pyrfume-data/chae_2019/pyrfumeData/' 

# mitral 33 data: 

filenames_mc33 = list(mitral_33.keys())
data_mc33 = list(mitral_33.values())

for i in range(len(filenames_mc33)):
    filename = filenames_mc33[i] + '_33.csv'
    data = pd.DataFrame(data = data_mc33[i])
    data.to_csv(header=False, index = False, path_or_buf = root + filename)

# mitral 55 data:
filenames_mc55 = list(mitral_55.keys())
data_mc55 = list(mitral_55.values())

for i in range(len(filenames_mc55)):
    filename = filenames_mc55[i] + '_55.csv'
    data = pd.DataFrame(data = data_mc55[i])
    data.to_csv(header=False, index = False, path_or_buf = root + filename)

# tufted data:
filenames_tufted = list(tufted.keys())
data_tufted = list(tufted.values())

for i in range(len(filenames_tufted)):
    filename = filenames_tufted[i] + '_55.csv'
    data = pd.DataFrame(data = data_tufted[i])
    data.to_csv(header=False, index = False, path_or_buf = root + filename)

