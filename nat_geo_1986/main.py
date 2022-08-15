#!/usr/bin/env python
# coding: utf-8

# # Processing workflow for the 1986 National Geographic Smell Survey data

# In[1]:


import pandas as pd
import pyrfume
from pyrfume.odorants import from_cids, get_cids


# In[2]:


df = pd.read_csv('NGS.csv', index_col=0).astype(int)  # Load the data
df.index.name = 'Subject'  # The index column is the subject number


# In[4]:


data_dict = pd.read_excel('Data_dictionary.xlsx', index_col=0)  # Load the data dictionary


# In[26]:


# Determine which integer value, if any, is used for no response (usually 0)
has_non_response_option = data_dict[data_dict['VALUES'].str.contains('No response') == True]
value_for_nan = has_non_response_option['VALUES'].apply(lambda x: x.split('=')[0]).astype(int)


# In[28]:


# Replace the value for no response with Python `None`)
# df = df.apply(lambda col: col.replace(value_for_nan.get(col.name, None), None))
df = df.apply(lambda col: col.replace(0, None)) # All "no reponse"'s are 0; bypasses regex error in above line


# In[29]:


# Odorant abbreviations used in the column names
odorant_abbreviations = {'AND': 'Androstenone',
                         'AA': 'Isoamyl acetate',
                         'AMY': 'Isoamyl acetate',
                         'GAL': 'Galaxolide',
                         'GALAX': 'Galaxolide',
                         'EUG': 'Eugenol', 
                         'MER': 'Mercaptans',
                         'MERCAP': 'Mercaptans',
                         'ROSE': 'Rose'}

# Question abbreviations used in the column names (see data dictionary for full question)
question_abbreviations = {'SMELL': 'Smell',
                          'QUAL': 'Quality',
                          'INT': 'Intensity',
                          'MEM': 'Memorable',
                          'EAT': 'Edible',
                          'WEAR': 'Wearable',
                          'DES': 'Descriptor'}
    
# List of unique odorant names
odorant_names = list(set(odorant_abbreviations.values()))


# In[30]:


# All (meta)-data not concerning the odorants themselves, i.e. information abou the subjects
metadata = df[[col for col in df if not any([col.startswith('%s_' % x) for x in odorant_abbreviations])]]
metadata.head()


# In[41]:


# Save this subject data
metadata.to_csv('subjects.csv')


# In[31]:


# All data concerning the odorants themselves
data = df[[col for col in df if any([col.startswith('%s_' % x) for x in odorant_abbreviations])]]

def f(s):
    """Convert e.g. 'AA_QUAL' into ('Amyl Acetate', 'Quality')"""
    odorant, question = s.split('_')
    return odorant_abbreviations[odorant], question_abbreviations[question]

# Turn column header into a multiindex with odorants names and questions as separate levels
data.columns = pd.MultiIndex.from_tuples(data.columns.map(f).tolist(), names=('Odorant', 'Question'))
data.head()


# In[32]:


# From methods.txt
# PEA added due to common knowledge that
# this is primary ingredient of IFF rose
molecule_names = ['5a-androst-16-en-3-one',
                  'isoamyl acetate',
                  'Galaxolide',
                  'eugenol',
                  'tert-butyl mercaptan',
                  'isopropyl mercaptan',
                  'n-propyl mercaptan',
                  'sec-butyl mercaptan',
                  'phenyl ethyl alcohol']


# In[33]:


# Get PubChem IDs for each odorant
names_to_cids = get_cids(molecule_names)


# In[34]:


# Generate information about molecules
cids = list(names_to_cids.values())
molecules = pd.DataFrame(from_cids(cids)).set_index('CID').sort_index()
molecules.head()


# In[35]:


names_to_cids


# In[26]:


# Save this molecule data
molecules.to_csv('molecules.csv')


# In[37]:


stimuli = pd.DataFrame(index=odorant_names, columns=[0]+cids)
# v/v * components ratios
stimuli.loc['Mercaptans'] = 0.04*pd.Series({6387: 0.76,
                                             6364: 0.18,
                                             7848: 0.04,
                                             10560: 0.02})
stimuli.loc['Androstenone'] = 0.001*pd.Series({6852393: 1})
stimuli.loc['Isoamyl acetate'] = 1*pd.Series({31276: 1})
stimuli.loc['Eugenol'] = 1*pd.Series({3314: 1})
stimuli.loc['Galaxolide'] = 0.425*pd.Series({91497: 1})

# Using common knowledge that IFF Rose is ~40% PEA
stimuli.loc['Rose'] = 0.8*pd.Series({6054: 0.4, 0: 0.6})
stimuli = stimuli.fillna(0)
stimuli.index.name = 'Simuli'
stimuli.to_csv('stimuli.csv')
stimuli.head()


# In[40]:


# Convert odorant names to PubChem IDs and pivot dataframe
data = data.stack('Question')
data = data.T.stack('Subject')#.reorder_levels([1, 0])
data.index.rename('Stimulus', level='Odorant', inplace=True)
data.head()


# In[42]:


# Save this behavioral response data
data.to_csv('behavior.csv')

