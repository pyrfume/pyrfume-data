#!/usr/bin/env python
# coding: utf-8

from pyrfume.odorants import from_cids, get_cids
import pandas as pd
import pickle

odorant_codes = []
experiments = [1, 2, 3]
for i in experiments:
    with open('data_odorset%d.pkl' % i, 'rb') as f:
        raw_data = pickle.load(f)
        assert len(raw_data['odor_names']) == len(set(raw_data['odor_1st']))
        assert len(raw_data['odor_names']) == len(set(raw_data['odor_2nd']))
        odorant_codes += raw_data['odor_names']

# From the paper
odorants = {'CinAld': 'Cinnamaldehyde',
            'EB': 'Ethylbutyrate',
            '22MBAcd': '2-Methylbutyric acid',
            '22DMBAcd': '2,2-Dimethylbutyric acid',
            'CPAcd': 'Cyclopentanecarboxylic Acid',
            '2Hep': '2-Heptanone',
            'IBAcd': 'Isobutyric acid',
            'IVAcd': 'Isovaleric acid',
            '3Hep': '3-Heptanone',
            'EB': 'Ethylbutyrate',
            'ValAcd': 'Valeric acid',
            '3MVAcd': '3-Methylvaleric acid',
            '33DMBAcd': '3,3-Dimethylbutyric acid',
            'Pinene': '(+)-alpha-Pinene', 
            'BzAld': 'Benzaldehyde',
            'IVAcd': 'Isovaleric acid',
            '4MVAcd': '4-Methylvaleric acid',
            'HexAcd': 'Hexanoic acid',
            'CinAld': 'Cinnamaldehyde',
            '3Hep': '3-Heptanone',
            'EB': 'Ethylbutyrate',
            '2MBAcd': '2-Methylbutyric acid',
            'MVT': 'Methylvalerate',
            'PPA': 'Propionic acid',
            'ButAcd': 'Butyric acid',
            '5M2H': '5-Methyl-2-Hexanone'}

cids = get_cids(list(odorants.values()))

molecules = pd.DataFrame(from_cids(list(cids.values()))).set_index('CID').sort_index()
molecules.to_csv('molecules.csv')
molecules.head()

all_data = []
for experiment_id in experiments:
    with open('data_odorset%d.pkl' % experiment_id, 'rb') as f:
        raw_data = pickle.load(f)
        
        lookup = pd.DataFrame(raw_data['odor_names'], columns=['odorant_code'])
        lookup['odorant_names'] = lookup['odorant_code'].apply(lambda x: list(map(odorants.get, x.split('_'))))
        lookup['cids'] = lookup['odorant_names'].apply(lambda x: list(map(cids.get, x)))
        lookup.index = lookup.index + 1
        
        del raw_data['odor_names']
        data = pd.DataFrame.from_records(raw_data).astype(int)
        data['cid1'] = data['odor_1st'].apply(lookup['cids'].xs).apply(sorted)
        data['cid2'] = data['odor_2nd'].apply(lookup['cids'].xs).apply(sorted)
        data['hash1'] = data[['cid1', 'conc_1st']].apply(lambda x: hash(str(tuple(x))), axis=1)
        data['hash2'] = data[['cid2', 'conc_2nd']].apply(lambda x: hash(str(tuple(x))), axis=1)
        data['mouse_id'] += 100*(experiment_id)
        all_data.append(data)
all_data = pd.concat(all_data).drop(['odor_1st', 'odor_2nd'], axis=1).fillna(0)

all_hashes = list(set(all_data['hash1']) | set(all_data['hash2']))
all_data['stimulus_1'] = all_data['hash1'].apply(all_hashes.index)+1
all_data['stimulus_2'] = all_data['hash2'].apply(all_hashes.index)+1

behavior = all_data.set_index(['stimulus_1', 'stimulus_2', 'delay', 'mouse_id'])[['response']]
behavior.to_csv('behavior.csv')
behavior.head()

z1 = all_data.set_index('stimulus_1')[['cid1', 'conc_1st']]
z1.columns = ['CIDs', 'conc']
z2 = all_data.set_index('stimulus_2')[['cid2', 'conc_2nd']]
z2.columns = ['CIDs', 'conc']
stimuli = pd.concat([z1, z2])
stimuli = stimuli[~stimuli.index.duplicated()].sort_index()
stimuli.index.name = 'stimulus'
stimuli.to_csv('stimuli.csv')
stimuli.head()

