#!/usr/bin/env python
# coding: utf-8

from functools import partial
import pandas as pd
from pyrfume import get_cids, from_cids

behavior_1 = pd.read_excel('DravnieksGrid.xlsx', sheet_name='App', index_col=0)

def fix(s):
    return s.replace("'", "").replace('Arnine', 'Amine')

behavior_1.index = behavior_1.index.map(fix)
behavior_1.columns = behavior_1.columns.map(fix)

behavior_2 = pd.read_excel('DravnieksGrid.xlsx', sheet_name='Use', index_col=0)#.set_index('CAS')
behavior_2.index = behavior_2.index.map(fix)
behavior_2.columns = behavior_2.columns.map(fix)

raw = pd.read_excel('Dravnieks_molecules.xlsx').set_index('Name')
raw['Conc'] = raw.index.map(lambda x: 'low' if 'low' in x.lower() else 'high')

raw.loc['Hexenal-trans1', 'CID'] = 5281168
raw.loc['Sandiff', 'CID'] = 103005
raw.loc['Tetraquinone', 'CID'] = 5424
raw.loc['PhenylEthanolhighconc', 'CID'] = 6054
raw.loc['PhenylEthanollowconc', 'CID'] = 6054
raw.loc['Diola', 'CID'] = 78484
#raw.loc['MethylAcetaldehydeDiAce', 'CID'] = 8503 # Speculative

raw.head()

behavior_1.index = behavior_1.index.map(lambda x: '_'.join([x, raw['Conc'].xs(x)]))
behavior_2.index = behavior_2.index.map(lambda x: '_'.join([x, raw['Conc'].xs(x)]))
behavior_1 = behavior_1[behavior_1.index.notnull()]
behavior_2 = behavior_2[behavior_2.index.notnull()]
#behavior_1.index = behavior_1.index.astype(int)
#behavior_2.index = behavior_2.index.astype(int)
behavior_1.index.name = 'Stimulus'
behavior_2.index.name = 'Stimulus'
behavior_1.to_csv('behavior_1.csv')
behavior_2.to_csv('behavior_2.csv')

behavior_2.head()

molecules = pd.DataFrame(from_cids(raw['CID'].dropna().tolist())).set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()]
molecules.head()

molecules.to_csv('molecules.csv')

stimuli = raw[['CAS', 'CID', 'Conc']].copy()
stimuli['Name'] = raw.index
stimuli.index = stimuli.apply(lambda x: '%s_%s' % (x['Name'], x['Conc']), axis=1)
stimuli.index.name = 'Stimulus'
stimuli.head()

assert not set(stimuli.index) - set(behavior_1.index)
assert not set(behavior_1.index) - set(stimuli.index)

stimuli.to_csv('stimuli.csv')
