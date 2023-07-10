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
# Dataset #1 comparisons
df1 = pd.read_csv('experiment1_comparisons.csv', header=0, index_col=0, names=['A', 'B', 'Similarity'])

df1_cids = pd.read_csv('experiment1_cids.csv', index_col=0)
df1_cids = df1_cids.applymap(lambda x: x.replace('[', '').replace(']', '').strip().replace(' ', ','))

df1.loc[:, ['A', 'B']] = df1.loc[:, ['A','B']].applymap(lambda x: df1_cids.loc[x]['Mixture Cids'])
df1.index.name = 'Pair'

print(df1.shape)
df1.head()

# +
# Dataset #2 comparisons
df2 = pd.read_csv('experiment2_comparisons.csv', header=0, index_col=0, names=['A', 'B', 'Similarity'])

df2_cids = pd.read_csv('experiment2_cids.csv', index_col=0)
df2_cids = df2_cids.applymap(lambda x: x.replace('[', '').replace(']', '').strip().replace(' ', ','))

df2.loc[:, ['A', 'B']] = df2.loc[:, ['A', 'B']].applymap(lambda x: df2_cids.loc[x]['Mixture Cids'])
df2.index = [df1.index[-1] + int(x) for x in df2.index]
df2.index.name = 'Pair'

print(df2.shape)
df2.head()

# +
# Dataset #3 comparisons
df3 = pd.read_csv('experiment3_comparisons.csv', header=0, index_col=0, names=['A', 'B', 'Similarity'])
df3.index = [df2.index[-1] + int(x) for x in df3.index]
df3.index.name = 'Pair'

print(df3.shape)
df3.head()
# -

# Combine all similarity data into single dataframe
df = pd.concat([df1, df2, df3])
df.rename(columns={'A': 'StimulusA', 'B': 'StimulusB'}, inplace=True)
df.head()

# +
# Physicochemical descriptors
physics = pd.read_csv('pcbi.1003184.s002.csv', index_col=0)
physics.index = physics.index.astype(int)

print(physics.shape)
physics.head()

# +
# Get list of all unique CIDs
cids1 = df1_cids['Mixture Cids'].apply(str.split, args=(',')).sum()
cids2 = df2_cids['Mixture Cids'].apply(str.split, args=(',')).sum()
cids3 = list(df3[['A', 'B']].values.ravel())
cids4 = physics.index.to_list()

cids = cids1 + cids2 + cids3 + cids4
cids = list(set(map(int, cids)))
# -

molecules_info = pyrfume.from_cids(cids)

# +
molecules = pd.DataFrame(molecules_info).set_index('CID').sort_index()

print(molecules.shape)
molecules.head()
# -

# Create dataframe for stimuli.csv; all simuli are CIDs
stimuli = pd.DataFrame(molecules.index.copy(), index=molecules.index.copy())
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to disk
molecules.to_csv('molecules.csv')
df.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')
physics.to_csv('physics.csv')
