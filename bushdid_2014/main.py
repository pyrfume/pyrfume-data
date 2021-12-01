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

# # Pyrfume data processing pipeline for Bushdid et al, 2014

import pandas as pd
from pyrfume.odorants import get_cids, from_cids

# ### Supplementary Table S1 containing information about the molecules

# Load and trim white space
s1 = pd.read_csv('tableS1.csv', encoding='latin1').iloc[:, :4].dropna()
s1['Odorant name'] = s1['Odorant name'].apply(str.strip)
s1.head()

# Get PubChem IDs from CAS registry numbers
cids = get_cids(s1['C.A.S.'])

# Which molecules could not be mapped to uniqye PubChem IDs?
s1[s1['C.A.S.'].isin(['87-44-5', '17369-59-4', '5655-61-8', '18172-67-3', '110-01-1'])]

# Add those CIDs manually
cids.update({
    '17369-49-4': 5373603, # celeriax is 3-propylidene phthalide
    '87-44-5': 5281515, # caryophyllene is beta-caryophyllene
    '5655-61-8': 6448, # Assumed racemic mixture for iso-bornyl acetate
    '18172-67-3': 440967, # laevo-beta-pinene is (-)-beta-pinene
    '110-01-1': 1127, # C.A.S. for thiolane was wrong in source file (Excel error)
})
cids = pd.Series(cids, name='CID').astype(int)

# Add the CIDs into the Table S1 dataframe
s1 = s1.join(cids, on='C.A.S.').set_index('CID')
s1 = s1.replace('110-01-1', '110-01-0') # Fix thiolane CID
s1.head()

# Get standardized information for these molecules
info = from_cids(cids.values)

# Join data and save.
molecules = pd.DataFrame(info).set_index('CID').join(s1)
molecules.to_csv('molecules.csv')

# ### Supplementary Table S2 containing information about the mixtures and the triangle tests that used them

s2 = pd.read_csv('tableS2.csv', encoding='latin1')
s2 = s2[s2['Test UID'].str.isnumeric() > 0]  # Keep only numeric Test UIDs
s2['Test UID'] = s2['Test UID'].astype(int)

# +
# Extract only the trial results for each Test UID and subject
behavior = pd.melt(s2, id_vars=['Test UID'], value_vars=['subject %d' % i for i in range(1, 27)], var_name='Subject', value_name='Correct')
behavior = behavior.replace({'subject %d' % i: i for i in range(1, 27)})
behavior = behavior.dropna().set_index(['Test UID', 'Subject']).sort_index()
behavior = behavior.replace({'wrong': False, 'right': True}).astype(bool)

# Save behavioral data
behavior.to_csv('behavior.csv')
behavior.head()

# +
# Get all the trial metadata (including mixture composition)
mixtures = s2.loc[:, :'subject 1'].iloc[:, :-1]
mixtures = mixtures.rename(columns={'Unnamed: %d' % i: 'Molecule %d' % (i-5) for i in range(6, 36)})

# Replace raw odorant names with CIDs
cid_map = {v.lower():k for k,v in molecules['Odorant name'].items()}
mol_columns = ["Molecule %d" % i for i in range(1, 31)]
mixtures[mol_columns] = mixtures[mol_columns].apply(lambda x: x.str.lower()).fillna(0)
# Manual fixes
mixtures = mixtures.replace(
    {'\xa01-isopropyl-4-methylbenzene': '1-isopropyl-4-methylbenzene',
     '4-methyl-3-penten-2-one': '4-methylpent-3-en-2-one'})
mixtures = mixtures.replace(cid_map)
mixtures[mol_columns] = mixtures[mol_columns].astype(int)

# Cleanup other names
mixtures['Stimulus dilution'] = mixtures['Stimulus dilution'].replace({'1/4': 0.25, '1/2': 0.5, 'not diluted': 1})
mixtures[['Components in mixtures', 'Components that differ']] = mixtures[['Components in mixtures', 'Components that differ']].astype(int)
mixtures = mixtures.set_index(['Test UID', 'Answer']).sort_index()

# Sanity check
for i, answer in enumerate(mixtures.index.get_level_values('Answer')):
    if i % 3 == 0:
        assert answer=='right'
    if i % 3 != 0:
        assert answer=='wrong'
mixtures.head(9)
# -

# Save mixtures data
mixtures.to_csv('mixtures.csv')


