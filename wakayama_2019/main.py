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

df = pd.read_csv('wakayama-intensity.txt', sep='\t')

cas_cids = get_cids(df['CAS'], kind='name')

manual_cids = {'80449-58-7': 225700,
               '124899-75-8': 106729,
               '58985-02-7': 10353,  
               '113889-23-9': 44153588}
cas_cids.update(manual_cids)

df['CID'] = df['CAS'].apply(cas_cids.get)
assert all(df['CID']>0)

df = df.set_index('CID').sort_index()
df.head()

molecules = pd.DataFrame(from_cids(df.index)).set_index('CID')

molecules[['CAS', 'Original Name']] = df[['CAS', 'Name']]
molecules.to_csv('molecules.csv')

behavior = df[['I_max', 'C', 'D']]
behavior.to_csv('behavior.csv')
