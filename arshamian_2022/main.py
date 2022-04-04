#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pyrfume

# Mapping between PubChem IDs and molecule names
cids = {379: 'Octanoic acid', 1183: 'Vanillin', 3314: 'Eugenol', 6054: '2-Phenylethanol', 6549: 'Linalool',
        7762: 'Ethyl butyrate', 8077: 'Diethyl disulfide', 10430: 'Isovaleric acid', 18827: '1-Octen-3-ol', 
        32594: '2-Isobutyl-3-methoxypyrazine'}

name_to_cid = dict((v,k) for k,v in cids.items())

# Order of odorants in 'All data intensity.txt'
order = {1: '1-Octen-3-ol', 2: '2-Isobutyl-3-methoxypyrazine', 3: '2-Phenylethanol', 4: 'Diethyl disulfide',
         5: 'Ethyl butyrate', 6: 'Eugenol', 7: 'Isovaleric acid', 8: 'Linalool', 9: 'Octanoic acid', 10: 'Vanillin'}

# get info from cids
info_dict = pyrfume.from_cids(list(cids.keys()))

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# Universal pleasantness data -> behavior_1.csv
behav1 = pd.read_csv('data\\Universal Pleasantness.csv', index_col=0)
behav1.head()

# Create unique ParticipantID for each Participant/Group
participants = sorted(list(set(list(zip(behav1['Group'], behav1['Participant'])))))
to_pid = {tup: i+1 for i, tup in enumerate(participants)}

# map OdorName to CID, Participant/Group to ParticipantID, and re-index
behav1['CID'] = behav1.apply(lambda row: name_to_cid[row['OdorName']], axis=1)
behav1['ParticipantID'] = behav1.apply(lambda row: to_pid[(row['Group'], row['Participant'])], axis=1)
behav1 = behav1.set_index(['CID','ParticipantID']).sort_index()
behav1.head()

# Dataframe for subjects.csv
subj_dict = dict((v,k) for k,v in to_pid.items())
subjects = pd.DataFrame.from_dict(subj_dict, orient='index', columns=['Group','Participant'])
subjects.index.names = ['ParticipantID']
subjects.head()

# All data intensity -> behavior_2.csv
behav2 = pd.read_csv('data\\All data intensity.txt', sep='\t')

# Replace Participant/Group with ParticipantID
behav2.rename(columns={'Participant number ': 'Participant','Language/Group':'Group'}, inplace=True)
behav2['ParticipantID'] = behav2.apply(lambda row: to_pid[(row['Group'], row['Participant'])], axis=1)

# Replace oder #'s with odor names
for col in behav2.columns[2:-1]:
    behav2[col].replace(order, inplace=True)

behav2 = behav2.set_index(['ParticipantID']).sort_index()
behav2.head()    

# Mainland lab ranking -> behavior_3.csv
behav3 = pd.read_csv('data\\UniversalPleasantnessMainlandLabRanking.csv')
behav3.head()

# convert rater names to raterID for anonymization
raters = list(set(behav3['Rater'].tolist()))
rater_to_id = {}
for i, rater in enumerate(raters):
    rater_to_id[rater] = i + 1

behav3['RaterID'] = behav3.apply(lambda row: rater_to_id[row['Rater']], axis=1)
behav3 = behav3.set_index(['CID','RaterID']).sort_index()
behav3.head()

# write to disk
molecules.to_csv('molecules.csv')
subjects.to_csv('subjects.csv')
behav1.drop(['OdorName', 'Group', 'Participant'], axis=1).to_csv('behavior_1.csv')
behav2.drop(['Group', 'Participant'], axis=1).to_csv('behavior_2.csv')
behav3.drop(['Rater', 'OdorKey', 'OdorName'], axis=1).to_csv('behavior_3.csv')