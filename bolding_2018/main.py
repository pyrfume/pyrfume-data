# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# + id="4oBNzBxS3dKI"
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as scio
from scipy.io import loadmat, whosmat
from scipy.sparse import csr_matrix, find
from typing import Any, Dict
import pyrfume


# + id="T-eBDls85c7Q"
# Utilities for extracting .mat and managing the resultant dicts

# Returns a matlab array as a nested dictionary
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

# Unpack a nested dictionary into a flat dict w composite key names
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

# Filter out dictionary keys that don't match a given suffix
def filter_dict_bykey(dictionary, keystring):
  filtered_dict = {key : value for key, value in dictionary.items() if key.find(keystring) >= 0}
  return filtered_dict

# Filter out dictionary key-value pairs for which there are no data
def filter_dict_bycontent(dictionary):
  filtered_dict = {key : value for key, value in dictionary.items() if len(value) > 0}
  return filtered_dict


# + id="P5gxukyh34Xb"
# Reading in data, which are contained in .efd files
path = 'untracked/'
files = [f for f in os.listdir(path) if f.endswith('.efd')]

# + id="_5VGQpY64HRe"
# Extract info about odor sets, recording location for each data file, 
# which are contained in the file 'ExperimetCatalog_Simul.txt'
experiment_info = pd.read_csv(path + 'ExperimentCatalog_Simul.txt', delimiter=' ')
experiment_info['filename'] = experiment_info['recording'].str.split('\\').str[-1]

all_data = []

# loop through all files in the directory
for f in files:
    (loading, region, filename) = experiment_info[
    experiment_info['filename'] == f.rstrip('.efd')][
        ['loading', 'region', 'filename']].values[0]

    data = load_matfile(path + f)

    # flatten the .mat array
    d = unpack_dict(data) # spike counts are as arrays here
    d = unpack_dict(d) # unpack spike counts across columns

    # filter out keys w/ irrelevant data and also keys lacking data 
    filtered = filter_dict_bykey(d, 'ValveSpikes_SpikesDuringOdor')
    filtered = filter_dict_bycontent(filtered)
    filtered_df = pd.DataFrame(filtered).T
    filtered_df.reset_index(inplace=True)

    # melt spike counts into long format
    melted_df = pd.melt(filtered_df, id_vars = 'index', value_vars=filtered_df.columns[1:], var_name='odor series', value_name='spike count')

    # populate df with columns for stim/subject identification later
    melted_df['recording'] = filename
    melted_df['region'] = region
    melted_df['loading'] = loading
    melted_df['odor idx'] = melted_df['index'].str.split('_').str[-2]
    melted_df['cell'] = melted_df['index'].str.split('_').str[-1]
    melted_df = melted_df[['recording', 'odor idx', 'cell', 'odor series', 'spike count', 'loading', 'region']]
    all_data.append(melted_df)

# + id="xdVAWQyk9rG0"
all_data = pd.concat(all_data)

# + id="VpZmjlRhAYMp"
# Odorant Panels used in the experiments. (A) is a concentration series, 
# and (C) is a set of 6 monomolecular odorants @ 0.3 v/v dilution. 
# See Bolding and Franks (2018) and crcns_pcs-1_data.pdf for a description of 
# indices and loadings. 

mono_odorants = [
    'ethyl butyrate', 
    '2-hexanone', 
    'isoamyl acetate', 
    'hexanal', 
    'ethyl tiglate', 
    'ethyl acetate'
    ]

concs = [0.03, 0.1, 0.3, 1]

# odorant panel A (2 conc. series and series of 6 odorants)
panelA = {
    '0' : ('mineral oil', 1),
    '1' : (mono_odorants[0], concs[0]),
    '2' : (mono_odorants[0], concs[1]),
    '3' : (mono_odorants[0], concs[2]),
    '4' : (mono_odorants[0], concs[3]),
    '5' : (np.nan, np.nan), # not defined for this paper
    '6' : (mono_odorants[1], concs[2]),
    '7' : (mono_odorants[2], concs[2]),
    '8' : (np.nan, np.nan), # 
    '9' : (mono_odorants[3], concs[0]),
    '10' : (mono_odorants[3], concs[1]),
    '11' : (mono_odorants[3], concs[2]),
    '12' : (mono_odorants[3], concs[3]),
    '13' : (np.nan, np.nan), # 
    '14' : (mono_odorants[4], concs[2]),
    '15' : (mono_odorants[5], concs[2])
}

# odorant panel C (6 monomolecular odorants @ 1 conc.)
panelC = {
    '0' : (np.nan, np.nan), # not defined for this paper
    '1' : (np.nan, np.nan), # 
    '2' : (np.nan, np.nan), # 
    '3' : (np.nan, np.nan), # 
    '4' : ('mineral oil', 1),
    '5' : (mono_odorants[3], concs[2]),
    '6' : (mono_odorants[1], concs[2]),
    '7' : (mono_odorants[2], concs[2]),
    '8' : (np.nan, np.nan), # 
    '9' : (mono_odorants[5], concs[2]),
    '10' : (mono_odorants[0], concs[2]),
    '11' : (mono_odorants[4], concs[2]),
    '12' : (np.nan, np.nan), # 
    '13' : (np.nan, np.nan), # 
    '14' : (np.nan, np.nan), # 
    '15' : (np.nan, np.nan), #    
}

# + id="wTWvADZvDUbC"
all_data['brain area'] = all_data['region'].map({'B' : 'bulb', 'P' : 'piriform'})

# Stimuli
stimuli = all_data[['odor idx', 'loading']]

# Subjects
subjects = all_data[['recording', 'cell', 'odor series', 'region']]

# + id="cH5VBlR1E6lH"
# Map the loadings and odor indices to (odor, conc.) tuples
stimuli.loc[stimuli['loading'] == 'A', 'odor'] = stimuli.loc[
    stimuli['loading'] == 'A', 'odor idx'].map(panelA)

stimuli.loc[stimuli['loading'] == 'C', 'odor'] = stimuli.loc[
    stimuli['loading'] == 'C', 'odor idx'].map(panelC)

stimuli[['odor name', 'concs']] = stimuli['odor'].apply(lambda x: pd.Series(x))

# + id="JNWZYtggu-8g"
# Generate stimulus IDs and index the 'stimuli' df on these
stimuli.drop('odor', axis=1, inplace=True)
stimuli = stimuli.drop_duplicates()
stimuli['stim ID'] = stimuli['odor idx'] + '_' + stimuli['loading']
stimuli = stimuli.set_index('stim ID')

# Generate subject IDs and index the 'subjects' df on these
subjects['sub ID'] = (subjects['recording'] + '_' 
                  + subjects['cell'] + '_'
                  + subjects['odor series'].astype(str) + '_' 
                  + subjects['region'])
subjects = subjects.drop_duplicates()
subjects.set_index('sub ID', inplace=True)

# + id="lf4fpeofy0iU"
data_merge_stim = pd.merge(all_data, stimuli.reset_index(), how='left', on=['odor idx', 'loading'])
data_merge_subjects = pd.merge(data_merge_stim, subjects.reset_index(), how='left', on=['recording', 'cell', 'odor series', 'region'])

# + colab={"base_uri": "https://localhost:8080/", "height": 260} id="OH9nRxXPz52Y" outputId="516010d0-caf6-4937-8f2b-b5afd109a18b"
# behavior
behavior = data_merge_subjects[['sub ID', 'stim ID', 'spike count']]
behavior.rename(columns={'sub ID': 'Subject', 'stim ID': 'Stimulus'}, inplace=True)
behavior = behavior.set_index('Stimulus')

behavior.head()

# + colab={"base_uri": "https://localhost:8080/", "height": 101, "referenced_widgets": ["ced5c31cd09b4e5ea20e7c4148bcb041", "d0d31870ae38488d9c32c8aab2b080de", "468cb10bbbda4574ada36d9330181c2b", "6d3c6b935c604ec9a600f2b0739b56aa", "2b3e3d6c3b30406bbf8ed9935bf5d113", "5f9a48e1f7b843b48b2e77991941156c", "dbf8660efcc943419c24cdd28db957d9", "0887e604482744b7af6056e651514e8a", "4d8d4789e65f497ab8e2be6d06d0eb28", "3c13a20b95ec4165b457ec9cfbe38a82", "c5d37b20c2f9413cbeb705b7225b2ed3", "6f105bdd73534a898be96837887f97f3", "fe086cba78c24c56b814b95d4241f654", "2d2d31cd0d364dd39711214bd7923438", "1861967b23af4a67bbb4175a34539b45", "6aa087faae4f411b8e91d2fbf40e4af8", "13838e3b7f1146728aa576220644b1a7", "86663ec82ac14324bf7d1e7e55bbd666", "9ad172a1b92f49fe9d169f4d87371da7", "a12e7903322d4d6faf261cf0ea8d1bee", "0b9395a8e3cc40b992a24a835cf85abe", "7baaa20a9b684e2f92dfcdd449208521"]} id="YP3sRwiMghP3" outputId="7588986d-7000-43a8-98c6-ed97a9914701"
odorants = pyrfume.get_cids(mono_odorants, kind='name')
cids = pd.Series(odorants)
molecules = pd.DataFrame(pyrfume.from_cids(cids.values)).set_index('CID').sort_index()

# + id="9nmzlINphTAc"
molecules.head()

# + id="YjyzoFxshe-W"
cids_merged = pd.merge(stimuli.reset_index(), molecules.reset_index(), how='left', left_on = 'odor name', right_on = 'name')
stimuli = cids_merged[['stim ID', 'odor idx', 'loading', 'odor name', 'CID', 'concs']]

# + id="5L4sRzNAlUqc"
stimuli.set_index('stim ID', inplace=True)
stimuli = stimuli[~stimuli['odor name'].isna()]

stimuli.head()

# + id="EaW8JbO3kTuk"
# Write the final .csv files
behavior.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')
subjects.to_csv('subjects.csv')
molecules.to_csv('molecules.csv')
