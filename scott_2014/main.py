#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
from pyrfume import get_cids
from pyrfume import from_cids

#Url + Convert First table into dataframe
url = 'https://www.jneurosci.org/content/34/6/2025/tab-figures-data'
tables = pd.read_html(url)
df = tables[0]

# Separate Dataframe into two tables based on 2-column format
df1 = df.iloc[:, [0,1,2,3]]
df2 = df.iloc[:, [4,5,6,7]]

#Rename columns of second dataframe and concatenate into 1 dataframe
df2 = df2.rename(columns={x: x.split('.')[0] for x in df2.columns})

behavior = pd.concat([df1, df2]).set_index('Number')

#Getting Molecules and CIDs
names = behavior['Odorant'].tolist()
cids = pd.Series(get_cids(names))

cids['Methyl pyridone'] = 12755
cids['Ethylacetoacetate'] = 8868
cids['Cyclohex anone']= 7967
cids['Cineole-1-8'] = 2758
cids['Methyl nonanate'] = 15606
cids['Cineole-1-4'] = 10106
cids['Methyloctanoate'] = 8091
cids['Cyclohexylacetate'] = 12146

info = from_cids(cids.values)

molecules = pd.DataFrame(info).set_index('CID').sort_index()
molecules = molecules[~molecules.index.duplicated()] # Remove duplicate rows
molecules.head()

behavior['CID'] = behavior['Odorant'].apply(cids.xs)
behavior = behavior.set_index('CID')
behavior.index.name = 'Stimulus'
behavior.drop('Odorant', axis=1, inplace=True)
behavior.head()

stimuli = pd.DataFrame(molecules.index, index=molecules.index)
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv', index=False)
stimuli.to_csv('stimuli.csv')

