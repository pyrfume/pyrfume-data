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

import pandas as pd
import pyrfume

# +
# Load raw data
data = pd.read_excel('mmc2.xlsx')
data['Name'] = data['Odorant'].apply(lambda x: x.split(' (')[0])

print(data.shape)
data.head()
# -

# Get CIDs
names = list(data['Name'].unique())
cids = pyrfume.get_cids(names, kind='name')

# +
# Manually add missing CIDs by manually searching PubChem and/or Google
manual_add = {
    'Bis(2-methyl-3-furyl)disulphide)': 526624,
    'a-Amylcinnamaldehyde dimethyl acetal': 6005821,
    'Chanel No 5': None, # Perfume
    'Olibanum Coeur MD': None, # Appears to be a substance
    'Decanoic_Acid': 2969,
    '2_coumaranone': 68382,
    'Axe': None # Assuming Axe is the brand and not 6,7-dimethoxy-~{N}-[(1~{R})-1-naphthalen-1-ylethyl]quinazolin-4-amine
}

for compound, cid in manual_add.items():
    cids[compound] = cid
# -

data['CID'] = data['Name'].map(cids)
data['Stimulus'] = data['Name'].apply(lambda x: x.replace(' ', '_'))
data.head()

# Get molecules info
molecules = pd.DataFrame(
    pyrfume.from_cids(
        [v for v in cids.values() if v]
    )
).set_index('CID')

molecules.sort_index(inplace=True)
molecules.head()

# +
# Behavior data
behavior = data[
    ['OR Name (Mouse OR Convention)', 'Minimum Activating Concentration [M]', 'log2(Fold Change)', 'PValue', 'FDR', 'Stimulus']
].copy()
behavior.rename(
    columns={'OR Name (Mouse OR Convention)': 'Subject', 'PValue': 'P Value', 'log2(Fold Change)': 'log2 (Fold Change)'},
    inplace=True
)
behavior = behavior.set_index(['Stimulus', 'Subject']).sort_index()

behavior.head()

# +
# Stimuli
stimuli = data.set_index('Stimulus')[['Odorant', 'CID']].copy()
stimuli.sort_index(inplace=True)

stimuli.head()

# +
# Subjects
subjects = data['OR Name (Mouse OR Convention)'].copy().to_frame()
subjects['Subject'] = subjects['OR Name (Mouse OR Convention)']
subjects = subjects.set_index('Subject').sort_index()

subjects.head()
# -

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')
subjects.to_csv('subjects.csv')
stimuli.to_csv('stimuli.csv')
