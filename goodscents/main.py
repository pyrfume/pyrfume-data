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

# # Preprocessing for Goodscents data

# ### Initial preprocessing previously done by Contrebande Labs
# ### Preprocessing here is to convert this to the Pyrfume standard format

from itertools import chain
import numpy as np
import pandas as pd
import pyrfume
from pyrfume.odorants import get_cids, from_cids, canonical_smiles, smiles_to_mol
from rdkit.Chem.Descriptors import MolWt
from tqdm.auto import tqdm

# Load the data previously processed by Google form the Leffingwell raw source file (not available here)
opl = pd.read_csv('data_rw_opl.csv')

# Load the data previously processed by Google form the Leffingwell raw source file (not available here)
odor = pd.read_csv('data_rw_odor.csv')

# Load the data previously processed by Google form the Leffingwell raw source file (not available here)
contents = pd.read_csv('data_rw_contents.csv')

# The CIDs given by Contrebande
cids_1 = opl['CID'].fillna(0).apply(lambda x: x[0] if isinstance(x, pd.Series) else x).fillna(0).astype(int)

# Export CAS for CIDs not found by Contrebande
opl[cids_1==0]['CAS Number'].to_csv('cas_out_2.csv', index=False, header=False)

# Now run cas_out_1.csv through the PubChem exchange service (https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi)
# Save the output and load it:
cas_to_cid = pd.read_csv('cas_to_cids_2.txt', sep='\t', index_col=0, header=None)[1]
cas_to_cid[''] = 0
# Map CAS numbers to CIDs
cids_2 = opl[cids_1==0]['CAS Number'].fillna('').apply(cas_to_cid.get, args=(0,))
# Pick the first CID if there are multiple hits
cids_2 = cids_2.apply(lambda x: x[0] if isinstance(x, pd.Series) else x).fillna(0).astype(int)

# Export names for things not yet found by Contrebande or CAS
opl[cids_1==0][cids_2==0]['Name'].to_csv('name_out_3.csv', index=False, header=False)

# Now run name_out_2.csv through the PubChem exchange service (https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi)
# Save the output and load it:
name_to_cid = pd.read_csv('name_to_cids_3.txt', sep='\t', index_col=0, header=None)[1]
name_to_cid[''] = 0
# Map names to CIDs
cids_3 = opl[cids_1==0][cids_2==0]['Name'].fillna('').apply(name_to_cid.get, args=(0,))
# Pick the first CID if there are multiple hits
cids_3 = cids_3.apply(lambda x: x[0] if isinstance(x, pd.Series) else x).fillna(0).astype(int)

# Found by each method
len(cids_1[cids_1>0]), len(cids_2[cids_2>0]), len(cids_3[cids_3>0])

opl = opl.copy() # Bypass ridiculous setting on a copy warning
# Apply CAS-based CIDs only to those CIDs not found by Contrebande
opl.loc[cids_1==0, 'CID'] = cids_2.copy()
# Apply name-based CIDs only to those CIDs not found by Contrebande or by CAS
opl.loc[cids_1==0, :].loc[cids_2==0, 'CID'] = cids_3.copy()
opl['CID'] = opl['CID'].astype(int)

# + tags=[]
cids = list(set(opl['CID']) - set([0]))
# Get standard information from PubChem for the CIDs that were found
molecules = from_cids(cids)
# Convert to a DataFrame
molecules = pd.DataFrame(molecules).set_index('CID').sort_index()
# -

molecules.to_csv('molecules.csv')
molecules.head()

tgsc_id = opl.set_index('CID')['TGSC ID']
identifiers = tgsc_id[tgsc_id.index>0].sort_index()
identifiers.to_csv('identifiers.csv')

delivery = odor.set_index('TGSC ID')[['Concentration %', 'Solvent']]
delivery.to_csv('delivery.csv')

z = odor.set_index('TGSC ID')['Tags'].fillna('[]').apply(lambda x: x.replace("'s'", "s'")) # Fix some apostrophe issues
behavior_sparse = z.apply(lambda x: set(eval(x)))
all_tags = set.union(*list(behavior_sparse.values))
behavior_sparse.to_csv('behavior_sparse.csv')

behavior = pd.DataFrame(index=behavior_sparse.index, columns=all_tags)
for tag in tqdm(all_tags):
    behavior[tag] = behavior_sparse.apply(lambda x: tag in x).astype(int)

behavior.to_csv('behavior.csv')
