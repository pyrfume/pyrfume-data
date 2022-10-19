#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
from rdkit.Chem import MolFromSmiles, rdMolDescriptors

# Load list of SMILES
smi = pd.read_csv('PrestwickChemLib.smi', header=None, names=['SMILES'])
smi.head()

# Fetch CIDs from SMILES
cids = pyrfume.get_cids(smi.SMILES.to_list(), kind='smiles')

# Use negative numbers for molecules with missing CIDs
i = -1
for k, v in cids.items():
    if not v:
        cids[k] = i
        i -= 1
        
smi['CID'] = smi.SMILES.map(cids)

# Get molecule info
molecules = pd.DataFrame(pyrfume.from_cids(smi.CID.to_list())).set_index('CID')

# Add in molecules that got assigned negative CIDs
for cid, sm in smi[smi.CID < 0][['CID', 'SMILES']].itertuples(index=False):
    mw = rdMolDescriptors.CalcExactMolWt(MolFromSmiles(sm))
    molecules.loc[cid] = [mw, sm, None, None]

molecules.sort_index(inplace=True)
molecules.head()

# Write to disk
molecules.to_csv('molecules.csv')

