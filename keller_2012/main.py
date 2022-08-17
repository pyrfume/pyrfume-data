#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pyrfume
from pyrfume.odorants import from_cids, get_cids

file_name = '12868_2012_2787_MOESM1_ESM.xlsx'

subjects = pd.read_excel(file_name, sheet_name='demographics', dtype={'income': str}).set_index('UID')
subjects.index.name = 'Subject'
subjects.head()

odorants = pd.read_excel(file_name, sheet_name='odours and sequence of stimuli')

replacements = {'propylen glycol': 'propylene glycol',
                '2-decenal': 'CC(=CCCC(=CC=O)C)C', # racemic 19801
                'citral': 'CC(=CCCC(=CC=O)C)C', # racemic 8843
                'isobornyl acetate': 'CC(=O)O[C@@H]1C[C@@H]2CC[C@@]1(C2(C)C)C', # 61061 matches isomeric SMILES on GoodScents
               }

odorants['name'] = odorants['odour'].replace(replacements)

# Add stimulus column; use vial # for intensity/pleasantness data; create stimulus ID for descriptors and threshold data
ip_data = odorants[odorants['measures obtained'] == 'intensity/pleasantness'].copy()
ip_data['vial #'] = ip_data['vial #'].astype(int)
ip_data['stimulus'] = 'IP' + ip_data['vial #'].astype(str)

d_data = odorants[odorants['measures obtained'] == 'descriptors'].copy()
d_data['tmp'] = range(1, len(d_data) + 1)
d_data['stimulus'] = 'D' + d_data.tmp.astype(str)

t_data = odorants[odorants['measures obtained'] == 'threshold'].copy()
t_data['tmp'] = range(1, len(t_data) + 1)
t_data['stimulus'] = 'T' + t_data.tmp.astype(str)

odorants.loc[ip_data.index, 'stimulus'] = ip_data['stimulus'].copy()
odorants.loc[d_data.index, 'stimulus'] = d_data['stimulus'].copy()
odorants.loc[t_data.index, 'stimulus'] = t_data['stimulus'].copy()

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')   
names_to_cids = get_cids(odorants['name'].values)

# Add CIDs to odorants
odorants['CID'] = odorants.name.map(names_to_cids).replace(0, np.nan)
odorants.head()

cids = [v for v in names_to_cids.values() if v]
molecules = pd.DataFrame(from_cids(cids)).set_index('CID').sort_index()

# Remove duplicates if any
molecules = molecules[~molecules.index.duplicated()]
molecules.head()

# Create dataframe for stimuli.csv
stimuli = odorants.drop('odour', axis=1).copy()
stimuli = stimuli.set_index('stimulus').sort_index()
stimuli = stimuli[['vial #', 'CID', 'CAS#', 'name', 'concentration', 'dilution', 'solvent', 'measures obtained', 'notes']]
stimuli.head()

# Threshold data -> behavior_1.csv
thresholds = pd.read_excel(file_name, sheet_name='thresholds').set_index('UID')
thresholds.columns = thresholds.columns.map(names_to_cids)
thresholds.reset_index(inplace=True)
thresholds = pd.melt(thresholds, id_vars='UID', var_name='CID', value_name='thresholds')
thresholds.rename(columns={'UID': 'subject'}, inplace=True)

cid_to_stimulus = dict(zip(odorants[odorants['measures obtained'] == 'threshold'].CID.astype(int),
                           odorants[odorants['measures obtained'] == 'threshold'].stimulus))

thresholds['stimulus'] = thresholds.CID.map(cid_to_stimulus)
thresholds = thresholds.drop('CID', axis=1).set_index(['stimulus', 'subject']).sort_index()

thresholds.head()

# Descriptor data -> behavior_2.csv
name_to_stimulus = dict(zip(odorants[odorants['measures obtained'] == 'descriptors'].name,
                           odorants[odorants['measures obtained'] == 'descriptors'].stimulus))

descriptors = pd.read_excel(file_name, sheet_name='descriptors', header=[0, 1], index_col=[0, 1])

# Drop two artifact columns from excel read
descriptors.drop(('pentadecalactone', 'lavender.1'), axis=1, inplace=True)
descriptors.drop(('androstenone', 'musty moldy.1'), axis=1, inplace=True)

descriptors.index.names = ['subject', 'rep']
descriptors.columns.names = ['mixture', 'descriptor']
mixture_names = descriptors.columns.levels[0]
mixture_stimulus = mixture_names.map(name_to_stimulus)
descriptors.columns = descriptors.columns.set_levels(mixture_stimulus, level=0)

# Convert from wide to long format
descriptors = pd.melt(descriptors, var_name=['stimulus', 'descriptor'], value_name='value', ignore_index=False)
descriptors.reset_index(inplace=True)
descriptors = descriptors.set_index(['stimulus', 'subject', 'rep']).sort_index()

descriptors.head()

# Intensity data -> behavior_3.csv
intensity = pd.read_excel(file_name, sheet_name='intensity and pleasantness', header=[1, 2], index_col=[0, 1]).iloc[:, :138]

def fix_indices(df):
    df.index.names = ['subject', 'rep']
    df.columns.names = ['conc', 'mixture']
    
    concs = df.columns.get_level_values(0)
    odours = df.columns.get_level_values(1)

    # Join duplicate mixtures
    odours = odours.map(lambda x: x.split('.')[0])
    
    # Fix mixture names
    replacements = {'propylen glycol': 'propylene glycol'}  
    odours = odours.map(lambda x: replacements.get(x, x))

    stim_ids = []
    for conc, odour in zip(concs, odours):
        conc = conc.replace('medium', 'high')  # Naming covention        
        pair = odorants[(odorants['measures obtained'] == 'intensity/pleasantness') & (odorants.odour == odour)].copy()
        pair = pair[['odour', 'concentration', 'dilution', 'solvent', 'stimulus']]
        pair = pair[~pair.duplicated(subset=['odour', 'concentration', 'dilution', 'solvent'])]
        
        # Two concentrations of everything except solvent
        assert (pair.shape[0] == 2) or (odour in ['paraffin oil', 'propylene glycol']), print(odour, pair.shape[0])
    
        if odour not in ['paraffin oil', 'propylene glycol']:
            stim = pair[pair['concentration'] == conc].stimulus.values[0]
        else:
            stim = pair.stimulus.values[0]
        stim_ids.append(stim)
        
    df.columns = pd.Index(stim_ids, name='stimulus')
    return df

intensity = fix_indices(intensity)

# Convert from wide to long format
intensity = pd.melt(intensity, var_name=['stimulus'], value_name='intensity', ignore_index=False)
intensity.reset_index(inplace=True)
intensity = intensity.set_index(['stimulus', 'subject', 'rep']).sort_index()

intensity.head()

# Pleasantness data -> behavior_4.csv
pleasantness = pd.read_excel(file_name, sheet_name='intensity and pleasantness', header=[1, 2], index_col=[0, 1]).iloc[:, 143:]
pleasantness = fix_indices(pleasantness)

# Convert from wide to long format
pleasantness = pd.melt(pleasantness, var_name=['stimulus'], value_name='pleasantness', ignore_index=False)
pleasantness.reset_index(inplace=True)
pleasantness = pleasantness.set_index(['stimulus', 'subject', 'rep']).sort_index()

pleasantness.head()

# Write to disk
molecules.to_csv('molecules.csv')
subjects.to_csv('subjects.csv')
stimuli.to_csv('stimuli.csv')
thresholds.to_csv('behavior_1.csv')
descriptors.to_csv('behavior_2.csv')
intensity.to_csv('behavior_3.csv')
pleasantness.to_csv('behavior_4.csv')

