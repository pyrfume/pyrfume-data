#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume

# Load list of odorants
odorants = pd.read_excel('Table1.xlsx', index_col=0)
odorants['Odorant Set'] = 'A'
odorants.iloc[17:34, odorants.columns.get_loc('Odorant Set')] = 'B'
odorants.iloc[34:, odorants.columns.get_loc('Odorant Set')] = 'C'
odorants.head()

# Get CIDs from names
names = odorants['Odorant name (abbreviation)'].dropna().str.split("(").str[0].str.strip().to_list()

cids = pyrfume.get_cids(list(set(names)), kind='name')

# Get info for molecules.csv
molecules = pd.DataFrame(pyrfume.from_cids(list(cids.values()))).set_index('CID').sort_index()

molecules.head()

# Create dataframe for stimuli.csv
stimuli = odorants.dropna(thresh=3, axis=0).drop('Literature', axis=1).copy().reset_index()

stimuli = stimuli.apply(lambda x: x.str.strip())
stimuli['Presenting set'] = stimuli['Presenting set'].fillna(method='ffill')
stimuli['Name'] = stimuli['Odorant name (abbreviation)'].str.split("(").str[0].str.strip()
stimuli['Abbreviation'] = stimuli['Odorant name (abbreviation)'].str.split("(").str[1].str.replace(')', '', regex=True)
stimuli['CID'] = stimuli.Name.map(cids)
stimuli.drop('Odorant name (abbreviation)', axis=1, inplace=True)

stimuli = stimuli.groupby(by=['Odorant Set', 'Abbreviation', 'Name', 'CID', 'Dilution factor (human)',
                   'Dilution factor (rat)', 'Odor descriptors'], as_index=False).agg(', '.join)

stimuli['Stimulus'] = stimuli[['Odorant Set', 'Abbreviation']].agg('_'.join, axis=1)
stimuli = stimuli.set_index('Stimulus').sort_index()
stimuli['Odor descriptors'] = stimuli['Odor descriptors'].str.lower().str.split(', ').apply(';'.join)

stimuli.head()

# Load similarity data
# For Set A:
#   LLGMNa = Using glomerular activity patterns evoked by isoamyl propionate diluted to 1/87 of saturated vapor over the neat
#   material
#   LLGMNb = Using glomerular activity patterns evoked by isoamyl propionate diluted to 1/10 of saturated vapor over the neat 
#   material
# For Set B:
#   LLGMNa = Using glomerular activity patterns evoked by butyl butyrate diluted to 1/10 of saturated vapor over the neat
#   material
#   LLGMNb = Using glomerular activity patterns evoked by butyl butyrate diluted to 1/190 of saturated vapor over the neat 
#   material
# For Set C:
#   LLGMNa = The only network used

data_dict = {}
for sheet in ['A', 'B', 'C']:
    data_dict[sheet] = pd.read_excel('Table2.xlsx', sheet_name='Set ' + sheet, skipfooter=2, index_col=0)
    if sheet == 'C':
        data_dict[sheet].rename(columns={'LLGMN': 'LLGMNa'}, inplace=True)

df = pd.concat(data_dict, axis=0).reset_index()
df.rename(columns={'level_0': 'Set', 'level_1': 'Abbr'}, inplace=True)
df.Abbr = df.Abbr.str.split("(").str[1].str.replace(')', '', regex=True)

df['Stimulus'] = df[['Set', 'Abbr']].agg('_'.join, axis=1).str.strip()
df = df.drop(['Set', 'Abbr'], axis=1).set_index('Stimulus').sort_index()
df = df.apply(lambda x: x.str.strip())
df.head()

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
df.to_csv('behavior.csv')

