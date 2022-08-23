#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
from rdkit.Chem import MolFromSmiles, rdMolDescriptors

# Load list of SMILES
raw = pd.read_csv('FragranceDB.smi', sep=' ', header=None, names=['SMILES', 'FrDB'])
raw.head()

# Fetch CIDs from SMILES
cids = pyrfume.get_cids(raw.SMILES.to_list(), kind='smiles')

# Use negative numbers for molecules with missing CIDs
i = -1
for k, v in cids.items():
    if not v:
        cids[k] = i
        i -= 1
        
raw['CID'] = raw.SMILES.map(cids)

# Get molecule info
molecules = pd.DataFrame(pyrfume.from_cids(raw.CID.to_list())).set_index('CID')

# Add in molecules that got assigned negative CIDs
for cid, smi in raw[raw.CID < 0][['CID', 'SMILES']].itertuples(index=False):
    mw = rdMolDescriptors.CalcExactMolWt(MolFromSmiles(smi))
    molecules.loc[cid] = [mw, smi, None, None]

molecules.sort_index(inplace=True)
molecules.head()

# Write to disk
molecules.to_csv('molecules.csv')

