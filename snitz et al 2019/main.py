#!/usr/bin/env python
# coding: utf-8

# Full reproducability is dependent on published data files to be downloaded and saved to ./published_data
# Note: odorants.csv had to be manually created by copy/pasting from Table 1 in online version of manuscript

import pandas as pd
import numpy as np
import re
import pyrfume

# info on odorants used in this study
# I had to manually add CID for IPM (the diluent)
odorants = pd.read_csv('published_data\\odorants.csv', index_col=0).dropna(subset=['CID'])

# mapping for odor codes to CID
codeToCID = dict(zip(odorants['Odor code'].str.split().str.join(''),odorants['CID'].astype('int')))

odorants.head()

# get info from cids
# drop rows with no CID since these represent mixtures
cids = odorants['CID'].dropna().astype('int').tolist()
info_dict = pyrfume.from_cids(cids)

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# data from SmellSpace participants
df1 = pd.read_excel('published_data\\bjz014_suppl_supplementary_datafile1.xlsx', engine='openpyxl')
df1 = df1.applymap(lambda s: s.upper() if type(s) == str else s) # set all odor codes to uppercase
df1.head()

# SmellSpace data -> behavior_1.csv
odorIdx1 = [df1.columns.get_loc(col) for col in df1 if col.startswith('Odor')]
colNames1 = ['CID', 'UID', 'OdorCode'] + [re.sub(r'\d+$', '', x) for x in list(df1)[odorIdx1[0]+1:odorIdx1[0]+14]]

data_dict1 = {}
n = 0
for row in df1.to_numpy():
    for i in odorIdx1:
        if row[i] in codeToCID: # skip codes indicating mixtures
            idx = [1] + [k for k in range(i,i+14)]
            data_dict1[n] = [codeToCID[row[i]]] + [row[k] for k in idx]
            n += 1

behav1 = pd.DataFrame.from_dict(data_dict1, orient='index', columns=colNames1).set_index(['CID', 'UID']).sort_index()
behav1.insert(1, 'SessionNumber', None) # Add blank column to be consitent with columns in control data
behav1.drop(['OdorCode'], axis=1).head()

# data from in-lab control group
df2 = pd.read_excel('published_data\\bjz014_suppl_supplementary_datafile2.xlsx', engine='openpyxl')

# Add user id's
df2.insert(0, 'UID', np.repeat(np.arange(len(df2)/2),2))
df2['UID'] = df2['UID'].astype(int)

df2.head()

# Control group data -> behavior.2.csv
odorIdx2 = [df2.columns.get_loc(col) for col in df2 if col.startswith('Odor')]
colNames2 = ['CID', 'UID', 'SessionNumber', 'OdorCode'] + [re.sub(r'\d+$', '', x) for x in list(df2)[odorIdx2[0]+1:odorIdx2[0]+14]]

data_dict2 = {}
n = 0
for row in df2.to_numpy():
    for i in odorIdx2:
        idx = [0,3] + [k for k in range(i,i+14)]
        data_dict2[n] = [codeToCID[row[i]]] + [row[k] for k in idx]
        n += 1

behav2 = pd.DataFrame.from_dict(data_dict2, orient='index', columns=colNames2).set_index(['CID', 'UID', 'SessionNumber']).sort_values(by=['CID', 'UID', 'SessionNumber'])
behav2.drop(['OdorCode'], axis=1).head()

# subjects.csv file to map UID's to age, gender where available
subj1 = df1[['UID', 'years Old']]
subj1.columns = ['UID', 'Age']

subj2 = df2[['UID', 'age', 'Gender']]
subj2.columns = ['UID', 'Age', 'Gender']

subjects = pd.concat([subj1, subj2.drop_duplicates()]).set_index('UID').sort_index()
subjects.head()

# write to disk
molecules.to_csv('molecules.csv')
behav1.drop(['OdorCode'], axis=1).to_csv('behavior_1.csv')
behav2.drop(['OdorCode'], axis=1).to_csv('behavior_2.csv')
subjects.to_csv('subjects.csv')