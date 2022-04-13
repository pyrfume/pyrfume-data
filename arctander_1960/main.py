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
from pyrfume.odorants import get_cids, from_cids, canonical_smiles

data = pd.read_csv('arctander.csv')

# +
data['SMILES'] = data['SMILES'].apply(lambda x: canonical_smiles(x) if isinstance(x, str) else None).str.strip()

def replace(s):
    if s:
        s = s.replace("’", "").replace('Δ', 'delta-')
    return s

data['ChemicalName'] = data['ChemicalName'].apply(lambda x: x if isinstance(x, str) else None).str.strip().apply(replace)
data.head()
# -

cid_col = 'new_CID'  # alternative is 'CID', which has fewer entries
data[cid_col] = data[cid_col].fillna(0)
needs_cids = data[data[cid_col]<=0].copy()

len(needs_cids)

cas_no_cids = needs_cids[needs_cids['CAS'].notnull()]['CAS'].str.strip()
cas_no_cids.to_csv('untracked/cas_no_cids.txt', header=None, index=False)
cas_no_cids  # Only one, aluminum sulfide solution (SID: 347706187)

# +
# Pointless if there are few CAS and all are known to not have CIDs (e.g. solutions not compounds)
# cids_from_cas = pd.read_csv('cids_from_cas.txt', sep='\t', index_col=0, header=None)[1]
# needs_cids.loc[cas_no_cids.index, cid_col] = cas_no_cids.apply(cids_from_cas.get, args=(0,)).fillna(0)
# (needs_cids[cid_col]==0).sum()
# -

names_no_cids = needs_cids[needs_cids['CAS'].isnull() & 
                           needs_cids['ChemicalName'].notnull()]['ChemicalName']
names_no_cids.to_csv('untracked/names_no_cids.txt', header=None, index=False)
len(names_no_cids)

# +
# Use PubChem Exchange to go from names_no_cids.txt to cids_from_names.txt
# -

cids_from_names = pd.read_csv('untracked/cids_from_names.txt', sep='\t', index_col=0, header=None)[1].dropna().astype(int)
needs_cids.loc[names_no_cids.index, cid_col] = names_no_cids.apply(cids_from_names.get, args=(0,)).fillna(0)
(needs_cids[cid_col]==0).sum()

smiles_no_cids = needs_cids[needs_cids['CAS'].isnull() & 
                            needs_cids['ChemicalName'].isnull() &
                            needs_cids['SMILES'].notnull()]['SMILES']
smiles_no_cids.to_csv('untracked/smiles_no_cids.txt', header=None, index=False)
len(smiles_no_cids)

# +
# Use PubChem Exchange to go from smiles_no_cids.txt to cids_from_smiles.txt

# +
# Pointless if there are no remaining SMILES without CIDs
#cids_from_smiles = pd.read_csv('cids_from_smiles.txt', sep='\t', index_col=0, header=None)[1].dropna().astype(int)
#needs_cids.loc[smiles_no_cids.index, cid_col] = smiles_no_cids.apply(cids_from_smiles.get, args=(0,)).fillna(0)
#(needs_cids[cid_col]==0).sum()
# -

# Set the CIDSs that were found
data.loc[needs_cids.index, cid_col] = needs_cids[cid_col].fillna(0)
cids = list(set(data[cid_col].values) - set([0]))

molecules = pd.DataFrame(from_cids(cids)).set_index('CID').sort_index()

identifiers = data[['ChemicalName', 'CAS', cid_col, 'ArctanderNum']].set_index('ArctanderNum')
identifiers = identifiers.sort_index()
identifiers.index.name = 'Stimulus'
identifiers.head()

molecules.to_csv('molecules.csv')
identifiers.to_csv('identifiers.csv')

sparse = data['ChastretteDetails'].fillna('').str.split()
sparse.index.name = 'Stimulus'
sparse.colums = 'Labels'
sparse.to_csv('behavior_1_sparse.csv')

all_labels = set.union(*sparse.apply(set).values)
dense = pd.DataFrame(index=sparse.index, columns=all_labels)
dense = dense.apply(lambda x: pd.Series({label: label in sparse[x.name] for label in x.index}),
                    axis=1)
dense = dense.astype(int).sort_index(axis=1)
dense.to_csv('behavior_1.csv')
dense.head()

desc = data[['Description']]
desc.to_csv('behavior_2.csv')
