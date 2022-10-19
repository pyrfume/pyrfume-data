#!/usr/bin/env python
# coding: utf-8

# Preprocessing for Goodscents data
# Initial preprocessing previously done by Contrebande Labs
# Preprocessing here is to convert this to the Pyrfume standard format

import numpy as np
import pandas as pd
import pyrfume
import os

# Load the data previously processed by Google form the Leffingwell raw source file (not available here)
opl = pd.read_csv('data_rw_opl.csv')
print(opl.shape[0], 'rows')
print(opl[opl['TGSC ID'].duplicated(keep=False)].shape[0], 'duplicate IDs')
opl.head()

opl2 = pd.read_csv('opl.csv')
print(opl2.shape[0], 'rows')
print(opl2[opl2['TGSC OPL ID'].duplicated(keep=False)].shape[0], 'duplicate IDs')
opl2.head()

# Load the data previously processed by Google form the Leffingwell raw source file (not available here)
odor = pd.read_csv('data_rw_odor.csv')
print(odor.shape[0], 'rows')
print(odor[odor['TGSC ID'].duplicated(keep=False)].shape[0], 'duplicate IDs')
odor.head()

odor = odor[['TGSC ID', 'Tags', 'Concentration %', 'Solvent']]
odor.sort_values(['TGSC ID', 'Concentration %', 'Solvent'], na_position='last', inplace=True)

# # Fix some apostrophe issues and convert to lists
odor.Tags = odor.Tags.fillna('[]').apply(lambda x: eval(x.replace("'s'", "s'")))

# Merge tags for duplicate TGSC ID's and remove duplicate descriptors
# Assumption here is that the duplicate ID's from data_rw_odor.csv are the result of listings for different sources
odor = odor.groupby('TGSC ID', dropna=False).agg({'Tags': 'sum', 'Concentration %': 'first', 'Solvent': 'first'}).reset_index()
odor.Tags = odor.Tags.apply(lambda x: list(set(x)))

print(odor.shape[0], 'rows left after merging duplicates')
odor.head()

df = odor.copy()

df['OPL ID'] = df['TGSC ID'].map(dict(zip(opl['TGSC ID'], opl['TGSC OPL ID'])))
df['CAS'] = df['TGSC ID'].map(dict(zip(opl['TGSC ID'], opl['CAS Number'])))
df['CID'] = df['TGSC ID'].map(dict(zip(opl['TGSC ID'], opl['CID']))).replace(np.nan, 0)
df['SMILES'] = df['OPL ID'].map(dict(zip(opl2['TGSC OPL ID'], opl2['SMILES'])))
df['Name'] = df['OPL ID'].map(dict(zip(opl2['TGSC OPL ID'], opl2['Common Name'])))
df['IUPAC Name'] = df['OPL ID'].map(dict(zip(opl2['TGSC OPL ID'], opl2['IUPAC Name'])))

# Use 'TGSC OPL ID' as index and stimulus ID since TGSC OPL ID is unique in opl.csv
df.drop_duplicates(subset=['OPL ID'], inplace=True) # Remove duplicates
df['OPL ID'].fillna(df.CAS, inplace=True)
df.set_index('OPL ID', inplace=True)
df.index.name = 'Stimulus'

print(df.shape[0], 'rows')
print(df[df.CID == 0].shape[0], 'missing CIDs')
df.head()

# Try fetching missing CIDs using SMILES
smi = df[df.CID == 0].SMILES.dropna()
cid_from_smi = pyrfume.get_cids(smi.to_list(), kind='smiles')

df.loc[smi.index, 'CID'] = df.loc[smi.index, 'SMILES'].map(cid_from_smi)
print(df[df.CID == 0].shape[0], 'missing CIDs')

# Now try finding remaining missing using CAS #'s
cas = df[df.CID == 0].CAS.dropna()
cid_from_cas = pyrfume.get_cids(cas.to_list())

df.loc[cas.index, 'CID'] = df.loc[cas.index, 'CAS'].map(cid_from_cas)
print(df[df.CID == 0].shape[0], 'missing CIDs')

# Now try finding remaining missing using names
# Searching by IUPAC names at this stage does not yeild any new CIDs
names = df[df.CID == 0]['Name'].dropna()
cid_from_name = pyrfume.get_cids(names.to_list(), kind='name')

df.loc[names.index, 'CID'] = df.loc[names.index, 'Name'].map(cid_from_name)
print(df[df.CID == 0].shape[0], 'missing CIDs')

# Manually add some of the remaining missing CIDs by searching PubChem 
# and drop remaining missing CIDS if also missing SMILES
# There are no missing CIDs after this step
df.loc['93905-03-4', 'CID'] = 33166
df.loc['27177-85-1', 'CID'] = 314293
df.loc['NF0133', 'CID'] = 11008539 
df = df[(df.CID != 0) & (~df.SMILES.isna())]

df.CID = df.CID.astype(int)
print(df.shape)
df.head()

# Get standard information from PubChem for the CIDs that were found
cids = list(set(df.CID.to_list()))

molecules = pd.DataFrame(pyrfume.from_cids(cids)).set_index('CID').sort_index()

molecules.head()

# Dataframe for stimuli.csv
stimuli = df[['TGSC ID', 'CID', 'Concentration %', 'Solvent']].copy()
stimuli.sort_index(inplace=True)
stimuli.head()

# Odor descriptor tags -> behavior.csv
behavior = df['Tags'].copy().to_frame().sort_index()
behavior['Descriptors'] = behavior.Tags.apply(lambda x: ';'.join(x))
behavior.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.drop('Tags', axis=1).to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')

