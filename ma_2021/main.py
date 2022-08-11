#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from pyrfume.odorants import get_cids, from_cids

odor_info = pd.read_excel('data_in_brief_V2.xlsx', sheet_name='odor information')
ratings = pd.read_excel('data_in_brief_V2.xlsx', sheet_name='individual data')

cas = odor_info['CAS.']
cids = get_cids(cas.values)

odor_info['CID'] = odor_info['CAS.'].apply(cids.get)
odor_info = odor_info.set_index('Odorant') # For easy map/apply lookup
odor_info = odor_info.rename(index={'ethyl 3-(methylsulfanyl)propan oate': 'ethyl 3-(methylsulfanyl)propanoate'})

molecules = pd.DataFrame(from_cids(list(cids.values()))).set_index('CID').sort_index()
molecules.head()

stimuli = odor_info[['CID', 'Cons.(mg/mL)', 'Solvent', 'Purity']]
stimuli = stimuli.rename(columns={'Cons.(mg/mL)': 'Concentration (mg/mL)'})
stimuli.index = stimuli.CID
stimuli.index.name = 'Stimulus'
stimuli.sort_index(inplace=True)
stimuli.head()

ratings['CID A'] = ratings['odor A'].apply(odor_info['CID'].xs)
ratings['CID B'] = ratings['odor B'].apply(odor_info['CID'].xs)
ratings['Repeat'] = ratings['Repeat'].fillna(1).astype(int)
ratings = ratings.rename(columns={'Sub.': 'Subject', 'Repeat': 'Rep', 'CID A': 'Stimulus A', 'CID B': 'Stimulus B'})
ratings.columns = ratings.columns.map(str.strip)
behavior = ratings.drop(['Trial', 'R-Trial', 'odor A', 'odor B'], axis=1)
behavior = behavior.set_index(['Stimulus A', 'Stimulus B', 'Subject', 'Rep']).sort_index()
behavior.head()

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
behavior.to_csv('behavior.csv')

