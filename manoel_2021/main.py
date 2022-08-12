#!/usr/bin/env python
# coding: utf-8

# # Preprocessing for Manoel et al, 2021

import pandas as pd
from pyrfume.odorants import from_cids, cids_to_smiles

# Load data about the odorants (names and PubChem IDs)
rosetta = pd.read_excel("SmilesInfo2.xlsx", engine='openpyxl')
threeletter_to_pubchem = dict(rosetta.iloc[:, -1:-3:-1].values)
threeletter_to_pubchem["MEN"] = threeletter_to_pubchem["+MEN"]
threeletter_to_pubchem = pd.Series(threeletter_to_pubchem, name='CID').to_frame()

# Load raw mouse data (individual mouse level)
raw = pd.read_csv("raw_behavioral_scores_mouse_73_odorants.csv", index_col=0, header=1).dropna()
raw.index.name = "odor"
raw = raw.join(threeletter_to_pubchem).sort_index()
print("Raw data has %d mouse/odorant combos" % raw.shape[0])
raw.head()

cids = raw['CID'].unique()
smiles = cids_to_smiles(cids)

molecules = pd.DataFrame(from_cids(cids)).set_index('CID').sort_index()
molecules.head()

behavior = raw.set_index('CID').sort_index().copy()
behavior.index.name = 'Stimulus'
behavior.head()

# Create dataframe for stimuli.csv; all simuli are CIDs
stimuli = pd.DataFrame(molecules.index, index=molecules.index)
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')

