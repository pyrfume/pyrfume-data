#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from pyrfume.odorants import get_cids, from_cids

raw = pd.read_csv('ifra-fragrance-ingredient-glossary---oct-2019.csv')
cas = raw['CAS number']

cids = get_cids(cas)

raw['CID'] = raw['CAS number'].apply(cids.get)
raw = raw[raw['CID'].notnull() & raw['CID']>0]

molecules = pd.DataFrame(from_cids(raw['CID'])).set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()] # Remove duplicates
molecules.head()

columns = ['CID', 'Primary descriptor', 'Descriptor 2', 'Descriptor 3']
behavior = raw[columns].set_index('CID')
behavior = behavior.rename(columns={'Primary descriptor': 'Descriptor 1'})
behavior.index.name = 'Stimulus'
behavior.head()

# Create dataframe for stimuli.csv; all simuli are CIDs
stimuli = pd.DataFrame(molecules.index, index=molecules.index)
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')

