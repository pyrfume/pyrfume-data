#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from pyrfume.odorants import from_cids

df1 = pd.read_csv('experiment1_comparisons.csv',
            header=0,index_col=0,names=['A','B','Similarity'])
df1_cids = pd.read_csv('experiment1_cids.csv', index_col=0)
df1_cids = df1_cids.applymap(lambda x:x.replace('[','').replace(']','').strip().replace(' ',','))
df1_cids
df1.loc[:, ['A','B']] = df1.loc[:, ['A','B']].applymap(lambda x:df1_cids.loc[x]['Mixture Cids'])
df1.head()

df2 = pd.read_csv('experiment2_comparisons.csv',
            header=0,index_col=0,names=['A','B','Similarity'])
df2_cids = pd.read_csv('experiment2_cids.csv', index_col=0)
df2_cids = df2_cids.applymap(lambda x:x.replace('[','').replace(']','').strip().replace(' ',','))
df2_cids
df2.loc[:, ['A','B']] = df2.loc[:, ['A','B']].applymap(lambda x:df2_cids.loc[x]['Mixture Cids'])
df2.head()

df3 = pd.read_csv('experiment3_comparisons.csv',
            header=0,index_col=0,names=['A','B','Similarity'])
df3.head()

df1.index.name = 'Pair'
df2.index.name = 'Pair'
df2.index = [200+int(x) for x in df2.index]
df3.index.name = 'Pair'
df3.index = [300+int(x) for x in df3.index]
df = pd.concat([df1, df2, df3])
df = df.rename(columns={'A': 'StimulusA', 'B': 'StimulusB'})
df.head()

cids1 = df1_cids['Mixture Cids'].apply(str.split, args=(',')).sum()
cids2 = df2_cids['Mixture Cids'].apply(str.split, args=(',')).sum()
cids3 = list(df3[['A', 'B']].values.ravel())

cids = cids1 + cids2 + cids3
cids = list(set(map(int, cids)))

molecules_info = from_cids(cids)

molecules = pd.DataFrame(molecules_info).set_index('CID').sort_index()
molecules.head()

# Create dataframe for stimuli.csv; all simuli are CIDs
stimuli = pd.DataFrame(molecules.index.copy(), index=molecules.index.copy())
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to disk
molecules.to_csv('molecules.csv')
df.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')

