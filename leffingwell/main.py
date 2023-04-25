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

# # Preprocessing for Leffingwell data

# ### Substantial preprocessing previously done by the team at Google; see leffingwell_readme.pdf
# ### Preprocessing here is to convert this to the Pyrfume standard format

from itertools import chain
import numpy as np
import pandas as pd
import pyrfume
from pyrfume.odorants import get_cids, from_cids, canonical_smiles, smiles_to_mol, hash_smiles
from rdkit.Chem.Descriptors import MolWt
from tqdm.auto import tqdm

# Load the data previously processed by Google form the Leffingwell raw source file (not available here)
raw = pd.read_csv('leffingwell_data.csv').set_index('smiles')

raw.head()

# Obtain the PubChem IDs -- ~100 of these ~3500 molecules cannot be found in PubChem
cids = pyrfume.get_cids(raw.index.to_list(), kind='smiles')

# Try to get missing CIDs from CAS
missing_cas = raw.loc[[k for k, v in cids.items() if not v]]['cas'].dropna()
cids2 = pyrfume.get_cids(missing_cas.to_list())

# +
# Add the PubChem ID column
raw['CID'] = raw.index.map(cids.get)
raw.loc[missing_cas.index, 'CID'] = raw.loc[missing_cas.index, 'cas'].map(cids2)

# Searching by chemical name doesn't find any new CIDs
# Use temporary negative numbers to fill in remaining missing CIDs
n = -1
for idx in raw[raw['CID'] == 0].index.to_list():
    raw.loc[idx, 'CID'] = n
    n -= 1
# -

# Canonicalize SMILES strings
raw.index = map(canonical_smiles, raw.index)

# Get standard information from PubChem for the CIDs that were found
all_cids = list(set([v for v in cids.values() if v] + [v for v in cids2.values() if v]))
info_dict = from_cids(all_cids)

# Convert to a DataFrame
info = pd.DataFrame(info_dict).set_index('CID').sort_index()
info.head()

# +
# Join the PubChem standard information with the original data
df = raw.join(info, on='CID', how='left')

# Those smiles associated with negative CID
empty_smiles = df[df['CID'] < 0].index

# Fill 'IsomericSMILES' column for molecules with no CID using original SMILES from index
df.loc[empty_smiles, 'IsomericSMILES'] = df.loc[empty_smiles].index

# Fill 'name' column for molecules with no CID using original `chemical_name`
df.loc[empty_smiles, 'name'] = df.loc[empty_smiles, 'chemical_name']

# No `IUPACName` will be computed for molecules with no CID

# Fill 'MolecularWeight' column for molecules with no CID using SMILES-based MW calculation
mols = smiles_to_mol(empty_smiles, max_attempts=1000)
mws = pd.Series({smiles: MolWt(mol) for smiles, mol in mols.items()})
df.loc[empty_smiles, 'MolecularWeight'] = mws[empty_smiles]
df['MolecularWeight'] = df['MolecularWeight'].astype(float)

# Replace negative integer CIDs with hash of SMILES
df.CID = df.apply(lambda row: hash_smiles(row['IsomericSMILES']) if row.CID < 0 else row.CID, axis=1)

# +
# Create dataframe that will generate molecules.csv and stimuli.csv
stimuli = df[list(info) + ['CID', 'cas']].set_index('CID').copy().sort_index()
stimuli.reset_index(inplace=True)

stimuli.index = stimuli.CID
stimuli.index.name = 'Stimulus'
stimuli = stimuli[~stimuli.index.duplicated()].sort_index() # Drop one duplicated row
stimuli.head()
# -

molecules = stimuli.drop('cas', axis=1).copy()
molecules = molecules.set_index('CID').sort_index()
molecules.head()

# Create the `behavior` dataframe containing the label data; all applicable labels are contained in the `Labels` column
behavior_sparse = df[['odor_data', 'odor_labels_filtered', 'CID']].set_index('CID').sort_index()
behavior_sparse.columns = ['Raw Labels', 'Labels']
behavior_sparse.index.name = 'Stimulus'
behavior_sparse = behavior_sparse[~behavior_sparse.index.duplicated()] # Drop one duplicated row
behavior_sparse.head()

# + tags=[]
# Create a dense version of the above; each label will have its own binary-valued column
behavior = behavior_sparse.copy()
# All the labels in the dataset
all_labels = set(chain.from_iterable(behavior['Labels'].squeeze().map(eval)))
for label in tqdm(all_labels):
    behavior[label] = behavior['Labels'].squeeze().apply(lambda x: label in eval(x)).astype(int)
behavior = behavior.drop(['Raw Labels', 'Labels'], axis=1).sort_index(axis=1)
# -

behavior.index.name = 'Stimulus'
behavior.head()

# Write files to disk
molecules.to_csv('molecules.csv')
behavior_sparse.to_csv('behavior_sparse.csv')
behavior.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')
