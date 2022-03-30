#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pyrfume

# load S1.csv data from tables/
s1 = pd.read_csv('tables\\S1.csv')

# only keep columns though "vapor.pressure.best.available"
idx = s1.columns.get_loc('vapor.pressure.best.available')
s1.drop(s1.iloc[:, idx+1:], inplace=True, axis=1)
s1.head()

# load S3.csv data from tables/ and reindex to SMILES
s3 = pd.read_csv('tables\\S3.csv').set_index('SMILES')
s3.head()

# get cids from SMILES
smiles = s1['SMILES'].tolist()
cids = pyrfume.get_cids(smiles)

# manually add the missing CIDs
man_add = {'CC(=O)OC[C@]1(O[C@H]2O[C@@H](COC(=O)C)[C@@H]([C@H]([C@@H]2OC(=O)C)OC(=O)C)OC(=O)C)O[C@H]([C@@H]([C@H]1OC(=O)C)OC(=O)C)COC(=O)C': -1,
           'O[C@@H]1[C@@H](O[C@H]2O[C@H](C(=O)O)[C@H]([C@H]([C@@H]2O)O)O)[C@@H](O[C@@H]([C@H]1O)C(=O)O)O[C@@H]1CC[C@@]2([C@@H](C1(C)C)CC[C@@]1([C@H]2C(=O)C=C2[C@@]1(C)CC[C@@]1([C@H]2C[C@](C)(CC1)C(=O)O)C)C)C': -2}

for name in man_add:
    cids[name] = man_add[name]

info_dict = pyrfume.from_cids(list(cids.values()))

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# function to replace SMILES with CID in s1 and s3 dataframes
def smile_to_cid(smile, cids):
    return cids[smile]

# reindex s1 dataframe for behavior_1.csv
s1['CID'] = s1.apply(lambda row: smile_to_cid(row['SMILES'], cids), axis=1)
s1.set_index('CID', inplace=True)
s1.sort_index(inplace=True)
s1.head()

# write to disk
molecules.to_csv('molecules.csv')
s1.drop(['SMILES'], axis=1).to_csv('behavior_1.csv')
s3.to_csv('behavior_2.csv')