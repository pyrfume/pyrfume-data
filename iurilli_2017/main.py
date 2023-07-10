# ---
# jupyter:
#   jupytext:
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

# + id="jNNzFzRoATgV"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as scio
from scipy.io import loadmat, whosmat
from scipy.sparse import csr_matrix, find
from typing import Any, Dict
import pyrfume

# + id="9IOBL9XgAaFA"
# Grab the data
files = [file for file in os.listdir() if '.mat' in file]


# + id="8zpdpzmNIvHm"
# Utility functions for handling and formatting .mat files
# load_matfile is from @esmitt on github

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
def filter_dict(dictionary, suffix):
    filtered_dict = {key : value for key, value in dictionary.items() if key.endswith(suffix)}
    return filtered_dict


# turn a sparse {0,1} matrix into a list of (row, column) pairs where
# matrix elements = 1
def unpack_sparse_matrix(sparse_matrix):
    indices = find(sparse_matrix)
    pairs = list(zip(indices[0], indices[1]))
    return pairs


# calculate delta firing rate w/ stimulation, averaged across trials
# function takes an array of (trial, spike time) pairs
def calculate_delta_FR(spike_list, t_baseline_start, t_stim, stim_duration):
    ms_in_sec = 1000
    spike_list = np.array(spike_list)
    if spike_list.shape[0] < 2 or spike_list.shape[1] < 2:
        return np.nan
    deltas = [] # delta firing rate across trials 
  # only grab trials where there was at least 1 spike (otherwise delta is undefined)
    trials = np.unique(spike_list[:,0])

  # calculate delta firing rate per trial:
    for t in trials:
    # find all spike times on a given trial 
        spiketimes = spike_list[:,1][np.where(spike_list[:,0] == t)]
    # sort spikes into those in pre vs. post-stim intervals
        rate_pre =  np.sum((spiketimes > t_baseline_start) & (spiketimes < t_stim)) / (t_stim - t_baseline_start)
        rate_post = np.sum((spiketimes > t_stim) & (spiketimes < t_stim + stim_duration)) / stim_duration
    # ignore if no spike ocurred in the pre or post intervals
        if not(rate_pre==0 and rate_post==0): 
            deltas.append(rate_post - rate_pre)
    mean_delta = ms_in_sec * np.array(deltas).mean()
    return mean_delta


# + id="MUSDpwvgDICM"
# Parse and extract data
# each .mat is turned into a flattened dict w/ composite keynames. Only spiking data are retained. 

conditions = [f.split('/')[-1].split('.')[0] for f in files] # brain structure + odorant panel
structures = [c.split('_')[0] for c in conditions] # brain structure only
all_data_dfs = []

for i in range(len(files)):
    # load .mat to a nested python dict
    dict_mat = load_matfile(files[i])
  
    # flatten the dict and convert sparse matrices into lists
    d = unpack_dict(dict_mat)
    d_filtered = filter_dict(d, '_spikeMatrix')

    # map all sparse spike matrices to (trial_num, spiketime) lists
    d_processed = {key: value for key, value in map(lambda item: (item[0], unpack_sparse_matrix(item[1])), d_filtered.items())}

    # reformat df to make Pyrfume-compatible
    df = pd.DataFrame(data=[d_processed.values()], columns=d_processed.keys()).T
    index_str = df.index.str.split('_').str
    p_str = df.index.str.split('_')
    df.rename(columns={0 : 'spikes'}, inplace=True)
    df['shank'] = index_str[3]
    df['cell'] = index_str[6]
    df['odor'] = index_str[8]
    df = df.drop(df[df['odor'] == 'spikeMatrix'].index)
    df['structure'] = structures[i]
    df['condition'] = conditions[i]
    df['ID'] = df['condition'] + '_' + df['shank'] + df['cell'] + df['odor']
    df = df[['ID', 'structure', 'shank', 'cell', 'odor', 'spikes']].set_index('ID')
    all_data_dfs.append(df)

# final concatenated df w/ all structures, conditions, cells, odors
data = pd.concat(all_data_dfs)

# calculate delta FR for all spike trains
t_baseline_start = 3000 # in ms 
t_stim = 4000 # in ms 
stim_duration = 2000 # in ms
data['delta_FR'] = data['spikes'].map(lambda x: calculate_delta_FR(x, t_baseline_start, t_stim, stim_duration))

# + colab={"base_uri": "https://localhost:8080/", "height": 236} id="s8qmDTcGcGYn" outputId="63f9cac2-715c-40e3-ae78-f4ce67c182a6"
data.head()

# + id="2CXh5fNPZAn2"
# track concentrations in two columns, for reasons noted below

data['odorant panel'] = data.index.str.split('_').str[1]
dilutions = [1, 0.1, 0.01, 0.001, 0.0001]

# concentration is handled in two ways, depending on which experiment series is being considered:
# 1)concentrations are 100 ppm for all odors that aren't part of the concentration series ('CS') experiments. 
# This is tracked in the column 'conc in ppm' (which is nan for experiments where it doesn't apply)
# 2) CS experiments are successive 10x dilutions across the odorant panels, for 5 concentrations. 
# This is tracked in the column 'conc dilution from 85 mM' (which is nan for experiments where it doesn't apply)

data['conc in ppm'] = (data['odorant panel'] == 'CS').replace({True : np.nan, False : 100})
data['conc dilution from 85 mM'] = data['odor']
data['conc dilution from 85 mM'] = data['conc dilution from 85 mM'].astype(int).apply(
    lambda x : dilutions[x % 5]) 

# mask-out the non-CS odorants w/ nans
data.loc[(data['odorant panel'] != 'CS'), 'conc dilution from 85 mM'] = np.nan

# + id="BGI6ywIaBHZE"
# Odor panels used in the experiment

# Odor panel 1 (N=15), as described in the methods section. 
# These are stimuili for the *_15.mat files
odorants_15 = {
    '0' : 'pentanal',
    '1' : 'hexanal',
    '2' : 'heptanal',
    '3' : 'heptanol',
    '4' : 'octanol',
    '5' : 'nonanol',
    '6' : 'phenetol',
    '7' : 'guaiacol',
    '8' : 'm-cresol',
    '9' : '2,4,5-trimethylthiazole',
    '10' : '4,5-dimethylthiazole',
    '11' : '4-methylthiazole',
    '12' : 'trimethylamine',
    '13' : 'isopentylamine',
    '14' : '2-phenyl-ethylamine' 
}

# Odor panel 2 (N=10), as described in the methods section.
# These are stimuli for the *_AA.mat files
odorants_AA = {
    '0' : '2,3,5-trimethyl-3-thiazoline',
    '1' : '2-methylbutyric acid',
    '2' : '2-propyltiethane',
    '3' : '3-mercapto-3-methuybutan-1-ol',
    '4' : 'isopentylamine',
    '5' : '2,3-butanedione',
  # Note: in the paper fig. legends, geraniol (6) comes after 2PET (7), but in the 
  # methods section, the order is reversed. We assumed the latter to be correct
  # until notified otherwise. 
    '6' : 'geraniol', 
    '7' : '2-phenylethanol',
    '8' : 'Peanut oil',
    '9' : 'Estrus female urine'
}

# Odor Panel 3 (N=15) described in the paper's methods section. 
# These are stimuli for the *_CS.mat files. The panel is a 5x concentration series for 
# 3 odorants (15 unique stimuli). Concentrations are tracked/disambiguated in the 
# column 'conc dilution from 85 mM', in the 'stimuli' dataframe, above.  
odorants_CS = {
    '0' : '2 phenylethanol', 
    '1' : '2 phenylethanol', 
    '2' : '2 phenylethanol',
    '3' : '2 phenylethanol',
    '4' : '2 phenylethanol',
    '5' : 'isoamylacetate', 
    '6' : 'isoamylacetate', 
    '7' : 'isoamylacetate', 
    '8' : 'isoamylacetate', 
    '9' : 'isoamylacetate', 
    '10' : '2,3,4-trimethyl-3-thiazoline', 
    '11' : '2,3,4-trimethyl-3-thiazoline', 
    '12' : '2,3,4-trimethyl-3-thiazoline', 
    '13' : '2,3,4-trimethyl-3-thiazoline', 
    '14' : '2,3,4-trimethyl-3-thiazoline', 
}

# Odor panel 4 (N=13), as described in the methods section. 
# These are stimuli for the *_natMix.mat files.    
odorants_natMix = {
    '0' : 'Sunflower Butter',
    '1' : 'Peanut Butter',
    '2' : 'Female Mouse Urine',
    '3' : 'Male Mouse Urine',
    '4' : 'Wolf Urine',
    '5' : 'Bobcat Urine',
    '6' : 'Lavander Flowers',
    '7' : 'Rose Oil',
    '8' : 'Coffee Beans',
    '9' : 'Mint Leaves',
    '10' : 'Hickory Chips',
    '11' : 'Clove Budes',
    '12' : 'Eucalyptus Oil'
}

# master list of all odorants across all panels
all_odorants = [list(odorants_15.values()), list(odorants_AA.values()), 
                list(odorants_CS.values()), list(odorants_natMix.values())]

all_odorants = np.array([element for sublist in all_odorants for element in sublist])
all_odorants = np.unique(all_odorants)

# find CIDs
odorants = pyrfume.get_cids(all_odorants, kind='name')

# de-orphan molecules that weren't found
odorants['3-mercapto-3-methuybutan-1-ol'] =  520682

# + colab={"base_uri": "https://localhost:8080/", "height": 236} id="jEfl7qEwpdUz" outputId="ec8cb445-2def-4f56-ae4a-065f57415fc2"
# stimuli
stimuli = data[['odorant panel', 'odor', 'conc in ppm', 'conc dilution from 85 mM']].drop_duplicates()
stimuli['stimuli'] = stimuli['odorant panel'] + '_' + stimuli['odor']
stimuli.set_index('stimuli', inplace=True)

# merge in the odor names
stimuli['odor name'] = stimuli['odor']

# ** Make sure you've created the stimulus dicts already ** 
stimuli.loc[(stimuli['odorant panel'] == '15'), 'odor name'] = stimuli.loc[
    (stimuli['odorant panel'] == '15'), 'odor name'].map(odorants_15)

stimuli.loc[(stimuli['odorant panel'] == 'AA'), 'odor name'] = stimuli.loc[
    (stimuli['odorant panel'] == 'AA'), 'odor name'].map(odorants_AA)

stimuli.loc[(stimuli['odorant panel'] == 'CS'), 'odor name'] = stimuli.loc[
    (stimuli['odorant panel'] == 'CS'), 'odor name'].map(odorants_CS)

stimuli.loc[(stimuli['odorant panel'] == 'natMix'), 'odor name'] = stimuli.loc[
    (stimuli['odorant panel'] == 'natMix'), 'odor name'].map(odorants_natMix)

stimuli['CID'] = stimuli['odor name'].map(odorants)
stimuli['CID'].replace({0 : np.nan}, inplace=True)
stimuli = stimuli[['odorant panel', 'odor', 'odor name', 'CID', 'conc in ppm', 'conc dilution from 85 mM']]
stimuli.head(5)

# + id="krO2RL1Rut9c"
# subjects
subjects = data[['structure', 'odorant panel', 'shank', 'cell']]
subjects['subjects'] = subjects['structure'] + '_' + subjects['odorant panel'] + '_' + subjects['shank'] + '_' + subjects['cell']
subjects.set_index('subjects', inplace=True)
subjects.drop_duplicates(inplace=True)
subjects.head(5)

# + id="yvEpGrYpPT0a"
# behavior
data_stim_merge = pd.merge(data, stimuli.reset_index(), how='left', on=['odorant panel', 'odor'])
data_subjects_merge = pd.merge(data_stim_merge, subjects.reset_index(), how='left', on=['structure', 'odorant panel', 'shank', 'cell'])
behavior = data_subjects_merge[['stimuli', 'subjects', 'spikes', 'delta_FR']]
behavior.set_index('stimuli', inplace=True)
behavior.head(10)

# + colab={"base_uri": "https://localhost:8080/", "height": 49, "referenced_widgets": ["f28d5df36f2242069c4e68cfc6d0f7b6", "2e82cdf9cb7647f78eef9c142ec43f14", "2406df02651d481e82aeeb79c756dd1c", "eb3db10ca4be4e6eb330cb30a44e8674", "b7cbe8491139487ebe5b5508db163762", "d42abffff07e45c087e95c83930c0174", "0ca06e563873479cb5ccdf98aeb73c78", "90a80350b37f436bb3946874a66c1abd", "ab04641a1b424e4cbfbe1e7a81fffb73", "046c7037903548b6af9aedfdb092df82", "b20e09e2608248b3aa18e41bbc086978"]} id="4RWxneuFVOQj" outputId="fd6d96bc-97b5-4d2d-e574-885c3d9a745c"
# molecules
cids = pd.Series(odorants)
molecules = pd.DataFrame(pyrfume.from_cids(cids.values))
molecules.sort_values(by='CID', inplace=True)
molecules.set_index('CID', inplace=True)

# + id="qmFzrDnSVvjs"
# write the behavior, subjects, stimuli, and molecules files
stimuli.to_csv('stimuli.csv')
subjects.to_csv('subjects.csv')
molecules.to_csv('molecules.csv')
# -

# behavior.csv file is too large for github. Using LFS-cache system.
actual_file = os.path.join('untracked','behavior.csv')
behavior.to_csv(actual_file)
hashed = pyrfume.lfs_hash(actual_file)
df_hashed_behavior = pd.DataFrame(['redirect', hashed])
df_hashed_behavior.to_csv('behavior.csv', header=False, index=False)
