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

# # Preprocessing for Manoel et al, 2021

import pandas as pd
from pyrfume.odorants import from_cids, cids_to_smiles

# +
# Load data about the odorants (names and PubChem IDs)
rosetta = pd.read_excel("SmilesInfo2.xlsx")
threeletter_to_pubchem = dict(rosetta.iloc[:, -1:-3:-1].values)
threeletter_to_pubchem["MEN"] = threeletter_to_pubchem["+MEN"]
threeletter_to_pubchem = pd.Series(threeletter_to_pubchem, name='CID').to_frame()

# Load raw mouse data (individual mouse level)
raw = pd.read_csv(
    "raw behavioral scores mouse 73 odorants.csv", index_col=0, header=1
).dropna()
raw.index.name = "odor"
raw = raw.join(threeletter_to_pubchem).sort_index()
print("Raw data has %d mouse/odorant combos" % raw.shape[0])
# -

cids = raw['CID'].unique()
smiles = cids_to_smiles(cids)

molecules = pd.DataFrame(from_cids(cids))
molecules.to_csv('molecules.csv')

behavior = raw.set_index('CID').join(pd.Series(smiles, name='SMILES')).set_index('SMILES')
behavior.to_csv('behavior.csv')
