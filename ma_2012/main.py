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
behav1 = df1[['CAS Number', 'Dorsal Response']]
behav1['CID'] = behav1['CAS Number'].map(cids)
behav1 = behav1.set_index(['CID']).sort_index()
behav1.head()

# DeltaF/F data in sd02.csv
expts = ['GIA0512', 'GIA0513', 'GIA0409', 'GIA12v2', 'GIA11-right',
        'GIA11-left', 'GIA8', 'GIA7', 'GIA6', 'GIA5', 'GIA4', 'GIA1']

# Mapping for odor abbreviations to CIDs
chemList = pd.read_excel('sd02.xls', sheet_name='ChemicalList', engine='xlrd')
abbr_to_cid = dict(zip(chemList['CAS Number'].map(cids), chemList['Abbreviation']))
chemList.head()

# Read xls sheets with DeltaF data
df2 = pd.read_excel('sd02.xls', sheet_name=expts, engine='xlrd')
for exp in df2:
    df2[exp].drop('Unnamed: 3', axis=1, inplace=True)
    df2[exp].set_index('Glom', inplace=True)

# Reshape into single dataframe -> behavior_2.csv
behav2 = pd.concat(df2, axis=1)

col_list=[]
for tup in behav2.columns:
    new_tup = (tup[0], tup[1].replace("'","").split('.')[0])
    col_list.append(new_tup)
behav2.columns = pd.MultiIndex.from_tuples(col_list)
behav2.head()

# Yellow highlights in sdo2.xls indicate 'movement artifact' and were not included in published analysis
# Create 'Artifact mask' of same shape as behavior such at 0 indicates 'artifact present' and 1 otherwise
wb = open_workbook('sd02.xls', formatting_info=True)

df_dict = {}
for exp in expts:
    sheet = wb.sheet_by_name(exp)
    header = sheet.row_values(0, start_colx=0, end_colx=None)
    header = [s.replace("'","").split('.')[0] for s in header]


    bgcol=np.zeros([sheet.nrows-1,sheet.ncols])
    for row in range(1, sheet.nrows): # skip 1st row as header
        for col in range(sheet.ncols):
            cif = sheet.cell_xf_index(row, col)
            iif = wb.xf_list[cif]
            cbg = iif.background.pattern_colour_index
            bgcol[row-1,col] = 0 if int(cbg) == 13 else 1
        
    # Add to dict as dataframe, match indexing to DeltaF dataframe
    df_dict[exp] = pd.DataFrame(bgcol, columns=header)
    df_dict[exp].index = range(1,len(df_dict[exp])+1)
    df_dict[exp].index.name = 'Glom'

# Remove two 'Glom' columns
for exp in df_dict:
    df_dict[exp].drop(['Glom',''], axis=1, inplace=True)

# Reshape into single dataframe with shape matching behavior_2.csv
df4 = pd.concat(df_dict, axis=1)
df4.head()

# Dataframe for odor concentration index for each experiment
# Concentration is presented as percentage of saturated vapor (S.V.)
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

df5 = pd.DataFrame.from_dict(conc, orient='columns')
df5.index.name = 'Concentration Index'
df5.head()

# Dataframe to map in-lab odor abbreviations to CID's
df6 = pd.DataFrame.from_dict(abbr_to_cid, orient='index', columns=['Odor Abbreviation'])
df6.index.name = 'CID'
df6 = df6.sort_index()
df6.head()

# write to disk
molecules.to_csv('molecules.csv')
behav1.drop('CAS Number', axis=1).to_csv('behavior_1.csv')
behav2.to_csv('behavior_2.csv')
df4.to_csv('artifact_mask.csv')
df5.to_csv('odor_conc.csv')
df6.to_csv('odor_abbr.csv')