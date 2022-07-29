#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
from xlrd import open_workbook
import pyrfume

# Dorsal response data in sd01.csv
df1 = pd.read_excel('sd01.xlsx', header=2, index_col=0, engine='openpyxl')
df1.head()

# Get CIDs from CAS
cas = df1['CAS Number'].tolist()
cids = pyrfume.get_cids(cas)

# Manually add missing CIDs
man_add = {'54830-99-8': 108630, '25340-17-4': 8657}

for k, v in man_add.items():
    cids[k] = v

# get info from cids
info_dict = pyrfume.from_cids(list(cids.values()))

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# Dataframe for dorsal response -> behavior_1.csv
d_resp = {'None': 0, '*': 1, '**': 2, '***': 3} # replace 0,*,**,*** notation with integers
behav1 = df1[['CAS Number', 'Dorsal Response']].copy()
behav1['CID'] = behav1['CAS Number'].map(cids) # convert CAS to CID
behav1['Dorsal Response'] = behav1['Dorsal Response'].map(d_resp) # convert respose notation to integers
behav1 = behav1.set_index(['CID']).sort_index()
behav1.head()

# DeltaF/F data is in sd02.csv, each experiment in it's own worksheet
# In-lab codes for each experiment:
expts = ['GIA0512', 'GIA0513', 'GIA0409', 'GIA12v2', 'GIA11-right',
        'GIA11-left', 'GIA8', 'GIA7', 'GIA6', 'GIA5', 'GIA4', 'GIA1']

# Mapping for odor abbreviations to CIDs
chemList = pd.read_excel('sd02.xls', sheet_name='ChemicalList', engine='xlrd')
chemList['CID'] = chemList['CAS Number'].map(cids)
abbr_to_cid = dict(zip(chemList['Abbreviation'], chemList['CID']))
chemList.head()

# Concentrations used in each experiment; units are percentage of saturated vapor (S.V.)
conc = {'GIA0512': {1: 2.5e-4, 2: 2.5e-3, 3: 2.5e-2},
        'GIA0513': {1: 2.5e-4, 2: 2.5e-3, 3: 2.5e-2},
        'GIA0409': {1: 2.5e-4, 2: 2.5e-3, 3: 2.5e-2},
        'GIA12v2': {1: 7.5e-5, 2: 2.5e-4, 3: 7.5e-4, 4: 2.5e-3, 5: 7.5e-3, 6: 2.5e-2},
        'GIA11-right': {1: 7.5e-5, 2: 2.5e-4, 3: 7.5e-4, 4: 2.5e-3, 5: 7.5e-3, 6: 2.5e-2},
        'GIA11-left': {1: 7.5e-5, 2: 2.5e-4, 3: 7.5e-4, 4: 2.5e-3, 5: 7.5e-3, 6: 2.5e-2},
        'GIA8': {1: 7.5e-6, 2: 2.5e-5, 3: 7.5e-5, 4: 2.5e-4, 5: 7.5e-4, 6: 2.5e-3, 7: 7.5e-3, 8: 2.5e-2}, 
        'GIA7': {1: 7.5e-6, 2: 2.5e-5, 3: 7.5e-5, 4: 2.5e-4, 5: 7.5e-4, 6: 2.5e-3, 7: 7.5e-3, 8: 2.5e-2},
        'GIA6': {1: 7.5e-6, 2: 2.5e-5, 3: 7.5e-5, 4: 2.5e-4, 5: 7.5e-4, 6: 2.5e-3, 7: 7.5e-3, 8: 2.5e-2}, 
        'GIA5': {1: 7.5e-6, 2: 2.5e-5, 3: 7.5e-5, 4: 2.5e-4, 5: 7.5e-4, 6: 2.5e-3, 7: 7.5e-3, 8: 2.5e-2}, 
        'GIA4': {1: 7.5e-6, 2: 2.5e-5, 3: 7.5e-5, 4: 2.5e-4, 5: 7.5e-4, 6: 2.5e-3, 7: 7.5e-3, 8: 2.5e-2}, 
        'GIA1': {1: 7.5e-6, 2: 2.5e-5, 3: 7.5e-5, 4: 2.5e-4, 5: 7.5e-4, 6: 2.5e-3, 7: 7.5e-3, 8: 2.5e-2} }

# Read xls sheets with DeltaF data
df2 = pd.read_excel('sd02.xls', sheet_name=expts, engine='xlrd')
for exp in df2:
    df2[exp].drop('Unnamed: 3', axis=1, inplace=True)
    df2[exp].columns = [c.replace("'","").split('.')[0] for c in df2[exp].columns]
    df2[exp] = pd.melt(df2[exp], ['Glom', 'X', 'Y'], var_name='OdorCode', value_name='DeltaF/F')
    df2[exp].set_index('Glom', inplace=True)

# Reshape into single long form dataframe
df3 = pd.concat(df2, axis=0)
df3.sort_index().head()

# Yellow highlights in sdo2.xls indicate 'movement artifact' and were not included in published analysis
# Create boolean 'artifact mask' of same shape as behavior_1 using 1: 'artifact present', 0: 'no artifact'
wb = open_workbook('sd02.xls', formatting_info=True)

artifact = {}
for exp in expts:
    sheet = wb.sheet_by_name(exp)
    header = sheet.row_values(0, start_colx=0, end_colx=None)
    header = [s.replace("'","").split('.')[0] for s in header]


    bgcol=np.zeros([sheet.nrows-1,sheet.ncols])
    for row in range(1, sheet.nrows): # skip 1st row as header
        for col in range(sheet.ncols):
            cif = sheet.cell_xf_index(row, col)
            iif = wb.xf_list[cif]
            cbg = iif.background.pattern_colour_index # this is cell background color
            bgcol[row-1,col] = 1 if int(cbg) == 13 else 0 # check if yellow, if so set to 1
        
    # Add to dict as dataframe, match indexing to DeltaF dataframe
    artifact[exp] = pd.DataFrame(bgcol, columns=header)
    artifact[exp].index = range(1,len(artifact[exp])+1)
    artifact[exp].index.name = 'Glom'

# Remove two 'Glom' columns
for exp in artifact:
    artifact[exp].drop(['Glom', '', 'X', 'Y'], axis=1, inplace=True)
    artifact[exp] = pd.melt(artifact[exp], var_name='OdorCode', value_name='Artifact', ignore_index=False)

# Reshape into single long form dataframe
df4 = pd.concat(artifact, axis=0)
df4.sort_index().head()

# Combine DeltaF/F and artifact mask for behavior_2.csv
behav2 = pd.concat([df3, df4['Artifact']], axis=1)
behav2.reset_index(inplace=True)
behav2.rename(columns={'level_0': 'Experiment'}, inplace=True)

# Add CIDs from in-lab odor codes
behav2['CID'] = behav2['OdorCode'].str[:-1].map(abbr_to_cid)
# 5 odor codes in deltaF/F data do not appear in chemical list, result in no CID, need to drop
behav2.dropna(axis=0, subset=['CID'], inplace=True)
behav2 = behav2.astype({'CID': int})

# Add concentration (% of S.V.)
def add_conc(row, conc):
    return conc[row['Experiment']][int(row['OdorCode'][-1])]

behav2['Conc (% of S.V.)'] = behav2.apply(lambda row: add_conc(row, conc), axis=1)
behav2.set_index(['CID', 'Experiment', 'Glom', 'Conc (% of S.V.)'], inplace=True)
behav2.sort_index(inplace=True)
behav2.head()

# Dataframe for identifiers.csv
ident = pd.DataFrame.from_dict({v: k for k, v in abbr_to_cid.items()}, orient='index', columns=['In-lab Abbreviation'])
ident.index.name = 'CID'
ident.sort_index(inplace=True)
ident.head()

# write to disk
molecules.to_csv('molecules.csv')
behav1.drop('CAS Number', axis=1).to_csv('behavior_1.csv')
behav2.drop('OdorCode', axis=1).to_csv('behavior_2.csv')
ident.to_csv('identifiers.csv')