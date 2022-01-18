# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Processing workflow for the 1986 National Geographic Smell Survey data

import pandas as pd
import pyrfume
from pyrfume.odorants import from_cids, get_cids

df = pd.read_csv('NGS.csv', index_col=0).astype(int)  # Load the data
df.index.name = 'Subject'  # The index column is the subject number

data_dict = pd.read_excel('Data dictionary.xlsx', index_col=0)  # Load the data dictionary

# Determine which integer value, if any, is used for no response (usually 0)
has_non_response_option = data_dict[data_dict['VALUES'].str.contains('No response') == True]
value_for_nan = has_non_response_option['VALUES'].apply(lambda x: x.split('=')[0]).astype(int)

# Replace the value for no response with Python `None`)
df = df.apply(lambda col: col.replace(value_for_nan.get(col.name, None), None))

# +
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
# -

# All (meta)-data not concerning the odorants themselves, i.e. information abou the subjects
metadata = df[[col for col in df if not any([col.startswith('%s_' % x) for x in odorant_abbreviations])]]
metadata.head()

# Save this subject data
metadata.to_csv('subjects.csv')

# +
# All data concerning the odorants themselves
data = df[[col for col in df if any([col.startswith('%s_' % x) for x in odorant_abbreviations])]]

def f(s):
    """Convert e.g. 'AA_QUAL' into ('Amyl Acetate', 'Quality')"""
    odorant, question = s.split('_')
    return odorant_abbreviations[odorant], question_abbreviations[question]

# Turn column header into a multiindex with odorants names and questions as separate levels
data.columns = pd.MultiIndex.from_tuples(data.columns.map(f).tolist(), names=('Odorant', 'Question'))
data.head()
# -

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

# Get PubChem IDs for each odorant
names_to_cids = get_cids(molecule_names)

# Generate information about molecules
cids = list(names_to_cids.values())
molecules = pd.DataFrame(from_cids(cids)).set_index('CID').sort_index()
molecules.head()

names_to_cids

# Save this molecule data
molecules.to_csv('molecules.csv')

mixtures = pd.DataFrame(index=odorant_names, columns=[0]+cids)
# v/v * components ratios
mixtures.loc['Mercaptans'] = 0.04*pd.Series({6387: 0.76,
                                             6364: 0.18,
                                             7848: 0.04,
                                             10560: 0.02})
mixtures.loc['Androstenone'] = 0.001*pd.Series({6852393: 1})
mixtures.loc['Isoamyl acetate'] = 1*pd.Series({31276: 1})
mixtures.loc['Eugenol'] = 1*pd.Series({3314: 1})
mixtures.loc['Galaxolide'] = 0.425*pd.Series({91497: 1})
# Using common knowledge that IFF Rose is ~40% PEA
mixtures.loc['Rose'] = 0.8*pd.Series({6054: 0.4, 0: 0.6})
mixtures = mixtures.fillna(0)
mixtures.to_csv('mixtures.csv')
mixtures.head()

# Convert odorant names to PubChem IDs and pivot dataframe
data = data.stack('Question')
data = data.T.stack('Subject')#.reorder_levels([1, 0])
data.head()

# Save this behavioral response data
data.to_csv('behavior.csv')
