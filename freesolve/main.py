#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
import pickle

# Load data
df = pd.read_csv('database.txt', sep=';', header=None, skiprows=3)

df.columns = ['compound id', 'SMILES', 'iupac name', 'experimental value (kcal/mol)', 'experimental uncertainty (kcal/mol)',
              'Mobley group calculated value (GAFF) (kcal/mol)', 'calculated uncertainty (kcal/mol)', 'experimental reference',
              'calculated reference', 'text notes']

df.SMILES = df.SMILES.str.strip()
df['iupac name'] = df['iupac name'].str.strip()

print(df.shape)
df.head()

# Fetch CIDs from SMILES
cids = pyrfume.get_cids(df.SMILES.to_list(), kind='smiles')

# Add to dataframe
df['CID'] = df.SMILES.map(cids)

# Look for missing CIDs by name
missing_by_name = df[df.CID == 0]

cids2 = pyrfume.get_cids(missing_by_name['iupac name'].to_list(), kind='name')

# Manually add missing
cids2['1,3-bis-(nitrooxy)butane'] = 101940534

df.loc[missing_by_name.index, 'CID'] = df.loc[missing_by_name.index, 'iupac name'].map(cids2)

# Get molecule info for molecules.csv
molecules = pd.DataFrame(pyrfume.from_cids(df.CID.to_list())).set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()] # Remove any duplciates
molecules.head()

# hydration free engergies -> physics.csv
physics = df[['CID', 'experimental value (kcal/mol)', 'experimental uncertainty (kcal/mol)',
              'Mobley group calculated value (GAFF) (kcal/mol)', 'calculated uncertainty (kcal/mol)']].copy()

physics = physics.set_index('CID').sort_index()
physics.head()

# Write to disk
molecules.to_csv('molecules.csv')
physics.to_csv('stimuli.csv')

