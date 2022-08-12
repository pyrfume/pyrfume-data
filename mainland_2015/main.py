#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
from pyrfume import get_cids
from pyrfume import from_cids

odorants = pd.read_csv('Odors.tsv', sep='\t')
odorants.head()

# Try to fetch CIDs first with SMILES
cids = pd.Series(get_cids(odorants['SMILES'], kind='SMILES'))

# Add CIDs to odorants dataframe
odorants['CID'] = odorants['SMILES'].apply(cids.xs).astype(int)

# Missing CIDs appear to be substances; assign negative numbers for CID
missing = odorants[odorants.CID == 0].CID.to_frame()
missing.CID = range(-1, -(missing.shape[0] + 1), -1)

odorants.loc[missing.index, 'CID'] = missing.CID
odorants = odorants.set_index('CID').sort_index()
odorants.head()

molecules = pd.DataFrame(from_cids(odorants.reset_index()['CID'])).set_index('CID').sort_index()

molecules.head()

# Odor #'s are unique, use as stimulus IDs
stimuli = odorants.drop('SMILES', axis=1).reset_index().copy()
stimuli.rename(columns={'Odor': 'Stimulus', 'CASRegistryNum': 'CASNo'}, inplace=True)
stimuli = stimuli.set_index('Stimulus').sort_index()
stimuli.head()

# Use OR # as Subject ID
subjects = pd.read_csv('Receptors.tsv', sep='\t')
subjects.index = subjects.OR
subjects.index.name = 'Subject'
subjects.head()

screening = pd.read_csv('PrimaryScreen.tsv', sep='\t')
screening.dropna(subset='Odor', axis=0, inplace=True) # Remove rows with no data
screening = screening.rename(columns={'Odor': 'Stimulus', 'OR': 'Subject'})
screening[['Stimulus', 'Subject']] = screening[['Stimulus', 'Subject']].astype(int)
screening = screening.set_index(['Stimulus', 'Subject']).sort_index()
screening.head()

data = pd.read_excel('sdata20152-s2.xls')
data = data.rename(columns={'Odor': 'Stimulus', 'OR': 'Subject'})
data = data.set_index(['Stimulus', 'Subject']).sort_index()

data.head()
|
# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
subjects.to_csv('subjects.csv')
screening.to_csv('behavior_1.csv')
data.to_csv('behavior_2.csv')

