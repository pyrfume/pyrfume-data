# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as scio
import os
from scipy.io import loadmat, whosmat
from scipy.sparse import csr_matrix, find
from typing import Any, Dict
import pyrfume

# returns a matlab array as a nested dictionary
def load_matfile(filename: str) -> Dict:
        def parse_mat(element: Any):
            # lists (1D cell arrays usually) or numpy arrays as well
            if element.__class__ == np.ndarray and element.dtype == np.object_ and len(element.shape) > 0:
                return [parse_mat(entry) for entry in element]

            # matlab struct
            if element.__class__ == scio.matlab.mio5_params.mat_struct:
                return {fn: parse_mat(getattr(element, fn)) for fn in element._fieldnames}

            # regular numeric matrix, or a scalar
            return element

        mat = scio.loadmat(filename, struct_as_record=False, squeeze_me=True)
        dict_output = dict()

        # not considering the '__header__', '__version__', '__globals__'
        for key, value in mat.items():
            if not key.startswith('_'):
                dict_output[key] = parse_mat(mat[key])
        return dict_output

# unpack a nested dictionary into a flat dict w composite key names
def unpack_dict(dictionary):
    unpacked_dict = {}

    for key, value in dictionary.items():
        if isinstance(value, dict):
            unpacked_value = unpack_dict(value)
            for nested_key, nested_value in unpacked_value.items():
                unpacked_dict[f"{key}_{nested_key}"] = nested_value
        elif isinstance(value, list):
            for index, item in enumerate(value):
                if isinstance(item, dict):
                    unpacked_value = unpack_dict(item)
                    for nested_key, nested_value in unpacked_value.items():
                        unpacked_dict[f"{key}_{index}_{nested_key}"] = nested_value
                else:
                    unpacked_dict[f"{key}_{index}"] = item
        else:
            unpacked_dict[key] = value

    return unpacked_dict

# filter out dictionary keys that don't match a given suffix
def filter_dict(dictionary, keystring):
    filtered_dict = {key : value for key, value in dictionary.items() if key.find(keystring) >= 0}
    return filtered_dict

path = 'untracked/'
files = [(path + file) for file in os.listdir(path)]

all_spike_data = []

for f in files:
    data = load_matfile(f)
    d = unpack_dict(data)
    df = filter_dict(d, '_ValveSpikes_SpikesDuringOdor') 

    # flatten remaining dicts of dicts
    odors_concs = unpack_dict(df) 
    odors_concs_cells = unpack_dict(odors_concs) # final flattened dict
    # clean up columns
    spike_data = pd.DataFrame(odors_concs_cells).T
    spike_data.reset_index(inplace=True)
    spike_data['index'] = spike_data['index'].str.lstrip('efd_ValveSpikes_SpikesDuringOdor_')
    spike_data = pd.melt(spike_data, id_vars=['index'], value_vars=list(range(15)), var_name='trial', value_name='num_spikes')
    spike_data[['odor', 'conc', 'cell']] = spike_data['index'].str.split('_', expand=True)
    spike_data['filename'] = f.split('.mat')[0].split('/')[-1]
    spike_data.set_index('index')
    spike_data = spike_data[['filename', 'odor', 'conc', 'cell', 'num_spikes']]
    all_spike_data.append(spike_data)

all_spike_data = pd.concat(all_spike_data)

stimuli = all_spike_data[['odor', 'conc']].drop_duplicates()

# stimuli used in the experiment
odor_map= {
    '0' : 'control (blank)',
    '1' : '2-hexanone',
    '2' : 'octanal',
    '3' : 'isoamyl acetate',
    '4' : 'ethyl butyrate',
    '5' : 'valeraldehyde',
    '6' : 'ethyl tiglate',
    '7' : 'acetophenone',
    '8' : 'gamma terpinene',
    '9' : 'ethyl acetate',
    '10' : 'methyl tiglate'
}

# get molecules from pubchem
odorants = pyrfume.get_cids(odor_map.values(), kind='name')

# map from odors to cids
odors_to_cids = {k : odorants[v] for k,v in odor_map.items()}
stimuli['CID'] = stimuli['odor'].map(odors_to_cids).replace(0, np.nan)
# generate stimulus IDs, and index on these
stimuli['stimulus'] = stimuli['odor'] + '_' + stimuli['conc']
#stimuli.set_index('stimulus', inplace=True)
stimuli['name'] = stimuli['odor'].map(odor_map)
stimuli.head(5)

# indexing recordings (different file names) as 'sessions', for labeling ease
filenames = all_spike_data['filename'].unique()
sessions = list(range(1, len(filenames) + 1))
sessions = [str(s) for s in sessions]
files_to_sessions = dict(zip(filenames, sessions))
subjects = all_spike_data[['filename', 'cell']]
subjects['experiment'] = subjects['filename'].map(files_to_sessions)
subjects['subject'] = subjects['experiment'] + '_' + subjects['cell']
subjects = subjects[['subject', 'experiment', 'cell', 'filename']]
subjects = subjects.drop_duplicates()
subjects.head(5)

spikes_stim_merge = pd.merge(all_spike_data, stimuli, how='left', left_on = ['odor', 'conc'], right_on=['odor', 'conc'])
spikes_stim_merge

spikes_subj_merge = pd.merge(spikes_stim_merge, subjects, how='left', left_on = ['filename', 'cell'], right_on=['filename', 'cell'])
behavior = spikes_subj_merge[['stimulus', 'subject', 'num_spikes']]
behavior = behavior.set_index('stimulus')
stimuli = stimuli.set_index('stimulus')
subjects = subjects.set_index('subject')

subjects['filename'] = subjects['filename'] + '.mat'

molecules = pd.Series(odorants)
molecules = pd.DataFrame(pyrfume.from_cids(molecules.values))
molecules = molecules.set_index('CID').sort_values(by='CID')
molecules

# write csv files
stimuli.to_csv('stimuli.csv')
subjects.to_csv('subjects.csv')
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')

