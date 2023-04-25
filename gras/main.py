# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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
import pyrfume
from pyrfume.odorants import hash_smiles
from rdkit.Chem import MolFromSmiles, rdMolDescriptors

# Load raw data
raw = pd.read_csv('GRAS.smi', header=None, names=['SMILES', 'CAS'], sep='\t')
raw.head()

# Get CIDs from SMILES
cids = pyrfume.get_cids(raw['SMILES'].to_list(), kind='SMILES', verbose=False)

# Try to find missing using CAS
raw['CID'] = raw.SMILES.map(cids)
cas_for_missing = raw[raw.CID == 0]
cids2 = pyrfume.get_cids(cas_for_missing.CAS.to_list())

# Manually add those still missing by searching PubChem; use temporary negative numbers for those not found
cids2['9005-67-8'] = 22833389
cids2['12/9/5550'] = -1
cids2['97593-31-2'] = -2
cids2['100085-39-0'] = -3

raw.loc[cas_for_missing.index, 'CID'] = raw.loc[cas_for_missing.index, 'CAS'].map(cids2)

molecules = pd.DataFrame(pyrfume.from_cids(raw.CID.to_list())).set_index('CID')

# +
# Add in molecules that got assigned negative CIDs
for cid, smi in raw[raw.CID < 0][['CID', 'SMILES']].itertuples(index=False):
    mw = rdMolDescriptors.CalcExactMolWt(MolFromSmiles(smi))
    molecules.loc[cid] = [mw, smi, None, None]

# Replace negative integer CID with hash of SMILES
molecules.index = molecules.apply(lambda row: hash_smiles(row['IsomericSMILES']) if row.name < 0 else row.name, axis=1)

# Remove any duplicates
molecules = molecules[~molecules.index.duplicated()].sort_index()

print(molecules.shape)
molecules.head()
# -

# Write to disk
molecules.to_csv('molecules.csv')
