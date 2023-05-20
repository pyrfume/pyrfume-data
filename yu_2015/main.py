# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

import pyrfume
import pandas as pd

# +
# Odorant response data from Supplemental Table 1
data = pd.read_excel('Table_S1.xlsx', header=0)
data['Domain'].ffill(inplace=True)

# Mapping of odor codes to names, from Table S1 legen
code_to_name = {
    'Oct': '1-octanol',
    'Car': '(−)-carvone',
    'Cou': 'coumarin',
    'Benzyl': 'benzyl acetate',
    'Allyl': 'allyl phenyl acetate'
}

data.rename(columns=code_to_name, inplace=True)

data.head()
# -

# Get CIDs
names = list(code_to_name.values())
cids = pyrfume.get_cids(names, kind='name')

# Manually add missing CID for (−)-carvone
cids['(−)-carvone'] = 439570

# Get molecule info
molecules = pd.DataFrame(pyrfume.from_cids(list(cids.values()))).set_index('CID').sort_index()
molecules.head()

# +
# Behavior
# Drop 'Total' since not a stimlus and can be regenerated easily by summing all 5 odorant responses
behavior = data.drop(['Domain', 'Total'], axis=1).rename(columns={'MOR256-3 ORs': 'Subject'}).copy()
# Reshape
behavior = pd.melt(behavior, id_vars='Subject', var_name='Stimulus', value_name='Normalized Response')
behavior.Stimulus = behavior.Stimulus.map(cids)
behavior = behavior.set_index(['Stimulus', 'Subject']).sort_index()

behavior.head()

# +
# Stimuli
stimuli = pd.DataFrame.from_dict(code_to_name, orient='index').reset_index()
stimuli.columns = ['In-lab Code', 'Name']
stimuli['CID'] = stimuli.Name.map(cids)
stimuli['Conc. (uM)'] = 300 # Response values in Table S1 are reported at 300 uM
stimuli['Stimulus'] = stimuli.CID
stimuli = stimuli.set_index('Stimulus').sort_index()

stimuli.head()

# +
# Subjects
subjects = data[['Domain', 'MOR256-3 ORs']].copy()
subjects['Subject'] = subjects['MOR256-3 ORs']
subjects = subjects.set_index('Subject')
subjects.columns = ['Domain', 'MOR256-3 ORs (WT or mutant)']

subjects.head()
# -

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
subjects.to_csv('subjects.csv')
behavior.to_csv('behavior.csv')
