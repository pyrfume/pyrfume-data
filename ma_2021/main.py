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

odor_info = pd.read_excel('data in brief V2.xlsx', sheet_name='odor information')
ratings = pd.read_excel('data in brief V2.xlsx', sheet_name='individual data')

cas = odor_info['CAS.']
cids = get_cids(cas.values)

odor_info['CID'] = odor_info['CAS.'].apply(cids.get)
odor_info = odor_info.set_index('Odorant') # For easy map/apply lookup
odor_info = odor_info.rename(index={'ethyl 3-(methylsulfanyl)propan oate': 'ethyl 3-(methylsulfanyl)propanoate'})

molecules = pd.DataFrame(from_cids(list(cids.values()))).set_index('CID').sort_index()
molecules.to_csv('molecules.csv')
molecules.head()

mixtures = odor_info.set_index('CID')[['Cons.(mg/mL)', 'Solvent', 'Purity']].sort_index()
mixtures = mixtures.rename(columns={'Cons.(mg/mL)': 'Concentration (mg/mL)'})
mixtures.to_csv('mixtures.csv')
mixtures.head()

ratings['CID A'] = ratings['odor A'].apply(odor_info['CID'].xs)
ratings['CID B'] = ratings['odor B'].apply(odor_info['CID'].xs)
ratings['Repeat'] = ratings['Repeat'].fillna(1).astype(int)
ratings = ratings.rename(columns={'Sub.': 'Subject',
                                  'Repeat': 'Rep'})
ratings.columns = ratings.columns.map(str.strip)
behavior = ratings.drop(['Trial', 'R-Trial', 'odor A', 'odor B'], axis=1)
behavior = behavior.set_index(['CID A', 'CID B', 'Subject', 'Rep']).sort_index()
behavior = behavior
behavior.to_csv('behavior.csv')
behavior.head()
