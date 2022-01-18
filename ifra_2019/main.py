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

raw = pd.read_csv('ifra-fragrance-ingredient-glossary---oct-2019.csv')
cas = raw['CAS number']

# + jupyter={"outputs_hidden": true} tags=[]
cids = get_cids(cas)
# -

raw['CID'] = raw['CAS number'].apply(cids.get)
raw = raw[raw['CID'].notnull() & raw['CID']>0]

molecules = pd.DataFrame(from_cids(raw['CID'])).set_index('CID')
molecules.to_csv('molecules.csv')

columns = ['CID', 'Primary descriptor', 'Descriptor 2', 'Descriptor 3']
behavior = raw[columns].set_index('CID')
behavior = behavior.rename(columns={'Primary descriptor': 'Descriptor 1'})
behavior.to_csv('behavior.csv')
