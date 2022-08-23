#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
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

# Manually add those still missing by searching PubChem; use negative numbers for those not found
cids2['9005-67-8'] = 22833389
cids2['12/9/5550'] = -1
cids2['97593-31-2'] = -2
cids2['100085-39-0'] = -3

raw.loc[cas_for_missing.index, 'CID'] = raw.loc[cas_for_missing.index, 'CAS'].map(cids2)

molecules = pd.DataFrame(pyrfume.from_cids(raw.CID.to_list())).set_index('CID')

# Add in molecules that got assigned negative CIDs
for cid, smi in raw[raw.CID < 0][['CID', 'SMILES']].itertuples(index=False):
    mw = rdMolDescriptors.CalcExactMolWt(Chem.MolFromSmiles(smi))
    molecules.loc[cid] = [mw, smi, None, None]

molecules.sort_index(inplace=True)
molecules.head()

# Write to disk
molecules.to_csv('molecules.csv')

