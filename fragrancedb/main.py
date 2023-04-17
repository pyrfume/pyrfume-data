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

# Load list of SMILES
raw = pd.read_csv('FragranceDB.smi', sep=' ', header=None, names=['SMILES', 'FrDB'])
raw.head()

# Fetch CIDs from SMILES
cids = pyrfume.get_cids(raw.SMILES.to_list(), kind='smiles')

# +
# Use temporary negative integers for molecules with missing CIDs
i = -1
for k, v in cids.items():
    if not v:
        cids[k] = i
        i -= 1
        
raw['CID'] = raw.SMILES.map(cids)
# -

# Get molecule info
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
