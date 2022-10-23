#!/usr/bin/env python
# coding: utf-8

# odorants_raw.csv was manually created by copy/pasting from Table 1 in online version of manuscript
# complete_odorants_csv.ipynb must then be run to clean odorants_raw.csv and add CID, CAS lists for mixtures

import pandas as pd
import numpy as np
import re
import pyrfume

# info on odorants used in this study
odorants = pd.read_csv('odorants.csv', index_col=0)      
odorants.head()

# get info from cids

# mapping for odor codes to CID
codeToCID = dict(zip(odorants.index, odorants['CID']))

# # remove mixtures
for key in [key for key in codeToCID if len(codeToCID[key].split(',')) > 1]: del codeToCID[key]                   

cids = list(codeToCID.values())
info_dict = pyrfume.from_cids(cids)

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules = molecules[~molecules.index.duplicated()]
molecules.head()

# data from SmellSpace participants
df1 = pd.read_excel('bjz014_suppl_supplementary_datafile1.xlsx', engine='openpyxl')
df1 = df1.applymap(lambda s: s.upper() if type(s) == str else s) # set all odor codes to uppercase
df1.rename(columns={'UID': 'Subject'}, inplace=True)
df1.head()

# SmellSpace data -> behavior_1.csv
odorIdx1 = [df1.columns.get_loc(col) for col in df1 if col.startswith('Odor')]
terms = [re.sub(r'\d+$', '', x) for x in list(df1)[odorIdx1[0]+1:odorIdx1[0]+14]]
colNames1 = ['Subject', 'Stimulus'] + terms

data_dict1 = {}
n = 0
for row in df1.to_numpy():
    for i in odorIdx1:
        idx = [1] + [k for k in range(i,i+14)]
        data_dict1[n] = [row[k] for k in idx]
        n += 1

behav1 = pd.DataFrame.from_dict(data_dict1, orient='index', columns=colNames1).set_index(['Stimulus', 'Subject']).sort_index()
behav1.head()

# data from in-lab control group
df2 = pd.read_excel('bjz014_suppl_supplementary_datafile2.xlsx', engine='openpyxl')

# Add user id's
df2.insert(0, 'Subject', np.repeat(np.arange(len(df2)/2),2))
df2['Subject'] = df2['Subject'].astype(int)
df2.head()

# Control group data -> behavior.2.csv
odorIdx2 = [df2.columns.get_loc(col) for col in df2 if col.startswith('Odor')]
colNames2 = ['CID', 'Subject', 'SessionNumber', 'Stimulus'] + terms

data_dict2 = {}
n = 0
for row in df2.to_numpy():
    for i in odorIdx2:
        idx = [0,3] + [k for k in range(i,i+14)]
        data_dict2[n] = [int(codeToCID[row[i]])] + [row[k] for k in idx]
        n += 1

behav2 = pd.DataFrame.from_dict(data_dict2, orient='index', columns=colNames2)    .set_index(['Stimulus', 'Subject', 'SessionNumber'])    .sort_values(by=['Stimulus', 'Subject', 'SessionNumber'])
behav2.drop(['CID'], axis=1, inplace=True)
behav2.head()

# subjects.csv file to map Subject's to age, gender where available
subj1 = df1[['Subject', 'years Old']]
subj1.columns = ['Subject', 'Age']

subj2 = df2[['Subject', 'age', 'Gender']]
subj2.columns = ['Subject', 'Age', 'Gender']

subjects = pd.concat([subj1, subj2.drop_duplicates()]).set_index('Subject').sort_index()
subjects.head()

# Dataframe for stimuli.csv - map in-lab Odor Abbreviations to CID
stimuli = pd.DataFrame.from_dict(codeToCID, orient='index', columns=['CID']).sort_index()
stimuli.index.name = 'Stimulus'
stimuli.head()

# write to disk
molecules.to_csv('molecules.csv')
behav1.to_csv('behavior_1.csv')
behav2.to_csv('behavior_2.csv')
subjects.to_csv('subjects.csv')
stimuli.to_csv('stimuli.csv')

