#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from pyrfume.odorants import get_cids, from_cids

df = pd.read_csv('wakayama-intensity.txt', sep='\t')

cas_cids = get_cids(df['CAS'])

manual_cids = {'80449-58-7': 225700,
               '124899-75-8': 106729,
               '58985-02-7': 10353,  
               '113889-23-9': 44153588}
cas_cids.update(manual_cids)

df['CID'] = df['CAS'].apply(cas_cids.get)
assert all(df['CID']>0)

df = df.set_index('CID').sort_index()
df.head()

molecules = pd.DataFrame(from_cids(df.index)).set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()] # Remove duplicates
molecules.head()

stimuli = df[['CAS', 'Name']].copy()
stimuli['CID'] = stimuli.index
stimuli.index.name = 'Stimulus'
stimuli.rename(columns={'Name': 'Original Name'}, inplace=True)
stimuli.head()

behavior = df[['I_max', 'C', 'D']].copy()
behavior.index.name = 'Stimulus'
behavior.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')

