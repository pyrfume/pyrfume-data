#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pyrfume

# load S1.csv data
s1 = pd.read_csv('S1.csv')

# only keep columns though "vapor.pressure.best.available"
idx = s1.columns.get_loc('vapor.pressure.best.available')
s1.drop(s1.iloc[:, idx+1:], inplace=True, axis=1)
s1.head()

# load S2.csv data
s2 = pd.read_csv('S2.csv').set_index('SMILES').sort_index()
s2.index.name = 'Stimulus'
s2.head()

# load S3.csv data
s3 = pd.read_csv('S3.csv').set_index('SMILES')
s3.drop('HAC', axis=1, inplace=True)
s3.index.name = 'Stimulus'
s3.head()

# get cids from SMILES
smiles = s1['SMILES'].tolist()
cids = pyrfume.get_cids(smiles, kind='smiles')

info_dict = pyrfume.from_cids(list(cids.values()))

# manually add the missing CIDs
man_add = [{'CID': -1,
            'MolecularWeight': 678.6,
            'IsomericSMILES': 'CC(=O)OC[C@]1(O[C@H]2O[C@@H](COC(=O)C)[C@@H]([C@H]([C@@H]2OC(=O)C)OC(=O)C)OC(=O)C)O[C@H]([C@@H]([C@H]1OC(=O)C)OC(=O)C)COC(=O)C',
            'IUPACName': None,
            'name': 'diastereomer of sucrose octaacetate'},
           {'CID': -2,
            'MolecularWeight': 822.4,
            'IsomericSMILES': 'O[C@@H]1[C@@H](O[C@H]2O[C@H](C(=O)O)[C@H]([C@H]([C@@H]2O)O)O)[C@@H](O[C@@H]([C@H]1O)C(=O)O)O[C@@H]1CC[C@@]2([C@@H](C1(C)C)CC[C@@]1([C@H]2C(=O)C=C2[C@@]1(C)CC[C@@]1([C@H]2C[C@](C)(CC1)C(=O)O)C)C)C',
            'IUPACName': None,
            'name': 'diastereomer of glycyron'}]

for d in man_add:
    cids[d['IsomericSMILES']] = d['CID']
    
info_dict += man_add

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# reindex s1 dataframe for behavior_1.csv
s1.set_index('SMILES', inplace=True)
s1.index.name = 'Stimulus'
s1.sort_index(inplace=True)
s1.head()

# Dataframe for stimuli.csv; use SMILES as stimulus ID
stimuli = pd.DataFrame.from_dict(cids, orient='index', columns=['CID']).sort_index()
stimuli.index.name = 'Stimulus'
stimuli.head()

# write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
s1.to_csv('behavior_1.csv')
s2.to_csv('behavior_2.csv')
s3.to_csv('physics.csv')

