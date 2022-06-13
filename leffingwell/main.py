# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3
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
from pyrfume.odorants import get_cids, from_cids, canonical_smiles, smiles_to_mol
from rdkit.Chem.Descriptors import MolWt
from tqdm.auto import tqdm

# Load the data previously processed by Google form the Leffingwell raw source file (not available here)
raw = pd.read_csv('leffingwell_data.csv').set_index('smiles')

# Obtain the PubChem IDs -- ~100 of these ~3500 molecules cannot be found in PubChem
cids = pyrfume.get_cids(raw.index, kind='smiles')

# Add the PubChem ID column
# 0 will be used to indicate molecules that could not be found in PubChem
raw['CID'] = raw.index.map(cids.get)

# Canonicalize SMILES strings
raw.index = map(canonical_smiles, raw.index)

# Get standard information from PubChem for the CIDs that were found
info_dict = from_cids(list(set(cids.values())))
# Convert to a DataFrame
info = pd.DataFrame(info_dict).set_index('CID')

# +
# Join the PubChem standard information with the original data
df = raw.join(info, on='CID', how='left')

# Those smiles associated with no CID
empty_smiles = df[df['CID']==0].index

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
# -

# Create the `molecules` dataframe containing basic information about the molecules
molecules = df[list(info) + ['CID', 'cas']].set_index('CID').copy().sort_index()
molecules.head()

# Create the `behavior` dataframe containing the label data; all applicable labels are contained in the `Labels` column
behavior_sparse = df[['IsomericSMILES', 'odor_data', 'odor_labels_filtered', 'CID']].set_index('CID').sort_index()
behavior_sparse.columns = ['IsomericSMILES', 'Raw Labels', 'Labels']

# + tags=[]
# Create a dense version of the above; each label will have its own binary-valued column
behavior = behavior_sparse.copy()
# All the labels in the dataset
all_labels = set(chain.from_iterable(behavior['Labels'].squeeze().map(eval)))
for label in tqdm(all_labels):
    behavior[label] = behavior['Labels'].squeeze().apply(lambda x: label in eval(x)).astype(int)
behavior = behavior.drop(['Raw Labels', 'Labels'], axis=1).sort_index(axis=1)
# -

# Write files to disk
molecules.to_csv('molecules.csv')
behavior_sparse.to_csv('behavior_sparse.csv')
behavior.to_csv('behavior.csv')
