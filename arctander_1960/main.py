#!/usr/bin/env python
# coding: utf-8
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import pandas as pd
from pyrfume.odorants import get_cids, from_cids, canonical_smiles

# Load raw data
data = pd.read_csv('arctander.csv')

data['SMILES'] = data['SMILES'].apply(lambda x: canonical_smiles(x) if isinstance(x, str) else None).str.strip()

def replace(s):
    if s:
        s = s.replace("’", "").replace('Δ', 'delta-')
    return s

data['ChemicalName'] = data['ChemicalName'].apply(lambda x: x if isinstance(x, str) else None).str.strip().apply(replace)
data.head()

cid_col = 'new_CID'  # alternative is 'CID', which has fewer entries
data[cid_col] = data[cid_col].fillna(0)
needs_cids = data[data[cid_col]<=0].copy()

len(needs_cids)

cas_no_cids = needs_cids[needs_cids['CAS'].notnull()]['CAS'].str.strip()
cas_no_cids  # Only one, ammonium sulfide solution (SID: 347706187)

# Pointless if there are few CAS and all are known to not have CIDs (e.g. solutions not compounds)

# Try to find missing CIDs form names
names_no_cids = needs_cids[(needs_cids[cid_col] == 0) & (needs_cids['ChemicalName'].notnull())]['ChemicalName']

len(names_no_cids)

cids_from_names = get_cids(names_no_cids.to_list(), kind='name', verbose=False)
print('%s found by name' % len([k for k, v in cids_from_names.items() if v]))

needs_cids.loc[names_no_cids.index, cid_col] = needs_cids.loc[names_no_cids.index, 'ChemicalName'].map(cids_from_names)
(needs_cids[cid_col]==0).sum()

# Try to find missing CIDs by SMILES
smiles_no_cids = needs_cids[(needs_cids[cid_col] == 0)
                            & (needs_cids['SMILES'].notnull())
                            & (needs_cids['SMILES'] != '')]['SMILES']

len(smiles_no_cids)

cids_from_smiles = get_cids(smiles_no_cids.to_list(), kind='SMILES', verbose=False)
print('%s found by SMILES' % len([k for k, v in cids_from_smiles.items() if v]))

needs_cids.loc[smiles_no_cids.index, cid_col] = needs_cids.loc[smiles_no_cids.index, 'SMILES'].map(cids_from_smiles)
(needs_cids[cid_col]==0).sum()

# Try to find missing CIDs by InChI Key
inchikey_no_cids = needs_cids[(needs_cids[cid_col] == 0) & needs_cids['InChiKey'].notnull()]['InChiKey']
len(inchikey_no_cids)

cids_from_inchikey = get_cids(inchikey_no_cids.to_list(), kind='inchikey', verbose=False)
print('%s found by InChI Key' % len([k for k, v in cids_from_inchikey.items() if v]))

# None found by InChI Key

# Set the CIDSs that were found
data.loc[needs_cids.index, cid_col] = needs_cids[cid_col].fillna(0)
cids = list(set(data[cid_col].values) - set([0]))

molecules = pd.DataFrame(from_cids(cids)).set_index('CID').sort_index()

print(molecules.shape)
molecules.head()

stimuli = data[['ChemicalName', 'CAS', cid_col, 'ArctanderNum']].set_index('ArctanderNum')
stimuli = stimuli.sort_index()
stimuli.index.name = 'Stimulus'
stimuli.head()

# Dict for standardizing labels
label_map = {
    'aldehidic': 'aldehydic',
    'almondy': 'almond',
    'aromatique': 'aromatic',
    'basalmic': 'balsamic',
    'camphor': 'camphoraceous',
    'citrusy': 'citrus',
    'ehtereal': 'ethereal',
    'foral': 'floral',
    'fruit': 'fruity',
    'gas': 'gassy',
    'herbaceous': 'herbal',
    'leathery': 'leather',
    'lillac': 'lilac',
    'mint': 'minty',
    'mushroomy': 'mushroom',
    'musk': 'musky',
    'must': 'musty',
    'nut': 'nutty',
    'peachy': 'peach',
    'peppery': 'pepper',
    'piney': 'pine',
    'root': 'rooty',
    'rosy': 'rose',
    'vanlilin': 'vanillin',
    'woney': 'winey',
    'wood': 'woody'
}

# +
sparse = data.set_index('ArctanderNum')[['ChastretteDetails']].rename_axis('Stimulus')
sparse.rename(columns=({'ChastretteDetails': 'Labels'}), inplace=True)
sparse['Labels'] = sparse['Labels'].fillna('').str.replace(',', '').str.split()

# Standardize lables and fix spelling errors
sparse['Labels'] = sparse['Labels'].apply(
    lambda x: [label_map[term] if term in label_map else term for term in x]
)

# Make sure no duplicate labels
sparse['Labels'] = sparse['Labels'].apply(lambda x: list(set(x)))

# To Series
sparse = sparse['Labels']
# -

all_labels = list(set.union(*sparse.apply(set).values))

dense = pd.DataFrame(index=sparse.index, columns=all_labels)
dense = dense.apply(lambda x: pd.Series({label: label in sparse[x.name] for label in x.index}),
                    axis=1)
dense = dense.astype(int).sort_index(axis=1)
print(dense.shape)
dense.head()

# +
# Convert lists to ;-separated strings
sparse = sparse.to_frame()
sparse['Labels'] = sparse['Labels'].apply(';'.join)

print(sparse.shape)
sparse.head()
# -

desc = data.set_index('ArctanderNum')[['Description']]
desc.index.name = 'Stimulus'
desc.head()

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
dense.to_csv('behavior_1.csv')
sparse.to_csv('behavior_1_sparse.csv')
desc.to_csv('behavior_2.csv')
