# -*- coding: utf-8 -*-
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

import pandas as pd
from pyrfume.odorants import get_cids, from_cids

# From the supplemental materials of:
# "An Algorithm for 353 Odor Detection Thresholds in Humans"
# Abraham et al, 2012 (Chemical Senses; 2011 online publication date)
df = pd.read_excel('ThresholdsAbraham2011.xls')

# Get CIDs for SMILES given in the original data file
smiles = df['SMILES'].dropna()
smiles_cids = get_cids(smiles, kind='SMILES')

# Use these CIDs where possible
df['CID'] = df['SMILES'].apply(smiles_cids.get, None)

# Replace typos and odd spellings with correct molecule names (whole names)
subs = {'lsobutylaldehyde': 'isobutyraldehyde',
        'n-Decylaldehyde': 'decanal',
        'Methyl sec.butyl ketone': '3-Methyl-2-pentanone',
        'Methyl tert.butyl ketone': 'Pinacolone',
        'a-Pinene': 'alpha-Pinene',
        'Butyl cellosolve  acetate': '2-Butoxyethanol acetate',
        '2-n-Buthoxyethanol': '2-butoxyethanol',
        '1-8 Cineole': 'eucalyptol',
        'n-Propy n-butyrate': 'Propyl butyrate',
        'D-3-carene': 'delta-3-carene'}
df['Substance'] = df['Substance'].replace(subs)

# Replace typos and odd spellings with correct molecule names (parts of names)
subs = {'.': '-',
        'alfa': 'alpha',
        'ß': 'beta',
        'mercaptane': 'mercaptan',
        'acryrale': 'acrylate',
        '- ': '-'}
for key, value in subs.items():
    df['Substance'] = df['Substance'].str.replace(key, value, regex=False)

# Get CIDS for molecule names that did not have SMILES (or whose SMILES could not be used)
names = df[df['CID'].isnull() | (df['CID']==0)]['Substance']
name_cids = get_cids(names, kind='name')

# Use these CIDs where CIDs could not be found previously
df.loc[names.index, 'CID'] = df['Substance'].apply(name_cids.get, None)

# Verify that a CID has been found for all molecules
assert all(df['CID']>0)

# Use the CID as the index and discard other identifiers from original dataset
df = df.set_index('CID').drop(['Substance', 'SMILES', 'MW'], axis=1)

# +
# There are some duplicate entries so average over duplicates and indicate where this has occurred
counts = df.groupby('CID')['Log (1/ODT)'].count()
behavior = df.groupby('CID').mean()
behavior['Duplicates'] = counts - 1

# Save this to the behavior file
behavior.to_csv('behavior.csv')

# +
# Get molecular data from PubChem (for consistency)
molecules = pd.DataFrame(from_cids(df.index)).set_index('CID')

# Save this to the molecules file
molecules.to_csv('molecules.csv')
