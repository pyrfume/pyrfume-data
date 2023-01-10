#!/usr/bin/env python
# coding: utf-8

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyrfume
from scipy.io import loadmat
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style('dark')
get_ipython().run_line_magic('matplotlib', 'inline')

# load raw data from .mat files
raw = loadmat('raw/odormappingdata_simple_05_2022.mat')

# .mat file contents (modified version of the description given at https://github.com/WachowiakLab/OdorMappingData):
# Variables: 
#  'allrespmats': Structure array with fields indexed in the following order:
#     (1) omp111L, (2) omp111R, (3) omp112L, (4) omp112R, (5) omp113L, (6) omp113R, (7) omp114L, (8) omp114R.
#     Variables included: 'respmat', 'Isparseness', 'numeffordorsperglom', 'roidynrange', 'normrespmat', 'dim', 'xnormrespmat',
#     'maxgloms', 'psparseness', 'numrespglomsperodor', 'rhocorrx', 'ROIPos', 'roimaxresp'
#     All odorants are in the same order for all response matrices (see ‘odornameslist’ variable). ROI positions 
#     indicate centroid of each ROI, after visual registration by aligning the midline and caudal sinus. Units are microns, 
#     reference (zero) is midline (for ‘Xpos’) and caudal sinus (for ‘Ypos’). Note that the Xposition for ROIs from the left OB
#     is negative relative to midline.
#  'allconcs_prep': Calculated final delivered concentration of each odorant, indexed in the same order as response matrices
#     and odornames list, for each mouse, in mols/L. Note that concentrations estimated from vapor pressure, calibrated air 
#     dilution, and liquid dilution, assuming ideal behavior.
#  'odornameslist': List of all 185 odorants (plus two vehicle controls), indexed in same order as odors in response spectrum 
#     matrices.
#  'respmatrix_ompXXXX_df': response matrices for each hemibulb

# Note that 'ROIPos' within the 'allrespmats' strucutre is in the Matlabe "Table" format which can't be loaded by
# scipy.io.loadmat; must export directly from Matlab to .csv.
# In Matlab used: writetable() for each ROIPos in the 'allrespmats' structure array

# Parse data
hemibulbs = ['111L', '111R', '112L', '112R', '113L', '113R', '114L', '114R']

# DeltaF response matrics; Row is RIO#, column is odor #.
resp_mat = {}
for ob in hemibulbs:
    resp_mat[ob] = raw[f'respmatrix_omp{ob}_dF']

# ROI coordinates
ROIPos = {}
for ob in hemibulbs:
    ROIPos[ob] = pd.read_csv('raw/' + ob + '_ROIPos.csv', index_col=0)
        
# Odor lists are identical, only need to grab once. First 185 are the odorants; last 2 are 'empty' (control) and the 
# solvent (triglyceride)
odor_list = [x[0][0] for x in raw['odornameslist']]
odorant_names = [' '.join(x.split(' ')[1:]) for x in odor_list]

# But only unique data relavtive to "fullinfo" .mat is the odorant concentration data
all_concs = pd.DataFrame(raw['allconcs_prep'])

# Heatmaps of response matrices
sns.set_style("ticks", {"xtick.bottom": True, "ytick.left": True})
fig, axes = plt.subplots(8, 1, figsize=(16, 32))
plt.xlabel('Odorant')

i = 0
for ob, im in resp_mat.items():
    if i == 7:
        ax = sns.heatmap(im, xticklabels=odorant_names, yticklabels=20, cbar=None, vmax=50, ax=axes[i])
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=5, ha='left')
    else:
        ax = sns.heatmap(im, xticklabels=[], yticklabels=20, cbar=None, vmax=50, ax=axes[i])
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, va='center')
    ax.set(ylabel=ob + ': Glom. (ROI#)')
    i += 1

# Get cids
cids = pyrfume.get_cids(odorant_names)

# Manually add the missing CIDs
man_add = {'(+-)-geosmin': 15559490,
           '5-methyl heptan-3-one oxime': 90787,
           '2-methoxy-3(5 or 6)-isopropylpyrazine': 71311442,
           'ethyl-2,5-dihydro-4-methylthiazole': 58165301,
           'empty': -1} # this is the control...

for name in man_add:
    cids[name] = man_add[name]

# Get molecule info from cids                   
info_dict = pyrfume.from_cids(list(cids.values()))

# Create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# Data for stimuli.csv
stimuli = pd.DataFrame(all_concs, copy=True)
stimuli.columns=['111', '112', '113', '114']
stimuli['CID'] = stimuli.index
stimuli['CID'] = stimuli['CID'].map({v: k for v, k in enumerate(odorant_names[:-2])}).map(cids)

# Convert to long form
stimuli = stimuli.melt(id_vars='CID', var_name='Mouse', value_name='Conc. (mols/L)')

# Create stimulus ID for each unique CID_concentration combo and index on this
stimuli['Stimulus'] = stimuli['CID'].astype(str) + '_' + stimuli['Conc. (mols/L)'].apply(lambda x: format(x, '.2e'))

# # add the blank ('empty') and solvent (triglyceride; CID = 5460048)
stimuli.loc[-1] = [-1, '111', None, 'empty'] # [CID, Mouse, Conc, Stimulus]
stimuli.loc[-2] = [-1, '112', None, 'empty']
stimuli.loc[-3] = [-1, '113', None, 'empty']
stimuli.loc[-4] = [-1, '114', None, 'empty']
stimuli.loc[-5] = [5460048, '111', None, '5460048_0']
stimuli.loc[-6] = [5460048, '112', None, '5460048_0']
stimuli.loc[-7] = [5460048, '113', None, '5460048_0']
stimuli.loc[-8] = [5460048, '114', None, '5460048_0']

# # Create dict to convert CID to Stimulus in behavior files
cid_to_stim = stimuli.set_index(['CID', 'Mouse']).drop('Conc. (mols/L)', axis=1).T.to_dict(orient='list')
cid_to_stim = {k: v[0] for k, v in cid_to_stim.items()}

stimuli.drop(['Mouse'], axis=1, inplace=True)
stimuli = stimuli.set_index(['Stimulus']).sort_index()
stimuli.drop_duplicates(inplace=True)
stimuli.loc['empty', 'CID'] = np.nan # Remove placeholder CID that was assigned to the control
stimuli.head()

# Data for subjects.csv
subj_dict = {}
for exp in resp_mat.keys():
    subj_dict[exp] = {'Mouse ID': exp[:-1], 'Hemisphere': exp[-1]}
    
# Add in ROI coordinates
subjects = pd.concat(ROIPos, axis=0).reset_index()
subjects.rename(columns={'level_0': 'Experiment', 'Row': 'ROI#'}, inplace=True)
subjects['ROI#'] = subjects['ROI#'].str.split(' ').str[1].astype(int)

subjects['Temp'] = subjects.Experiment.map(subj_dict)
subjects['Mouse ID'] = subjects.Temp.apply(lambda d: d['Mouse ID'])
subjects['Hemisphere'] = subjects.Temp.apply(lambda d: d['Hemisphere'])
subjects['Subject'] = subjects['Experiment'] + '_' + subjects['ROI#'].astype(str).str.zfill(3)

subjects = subjects.set_index('Subject').sort_index()
subjects = subjects[['Experiment', 'Mouse ID', 'Hemisphere', 'ROI#', 'Xpos', 'Ypos']]
subjects.head()

# Create dataframe from response matrices -> behavior.csv
resp_mat_dict = {}
for key, mat in resp_mat.items():
    tmp = pd.DataFrame(mat)
    tmp.columns = [cids[odorant_names[i]] for i in tmp.columns]
    tmp['Experiment'] = key
    tmp.reset_index(inplace=True)
    tmp['index'] = tmp['index'] + 1
    tmp = tmp.melt(id_vars=['index', 'Experiment'], var_name='CID', value_name='DeltaF', ignore_index=False)
    tmp.rename(columns={'index': "ROI#"}, inplace=True)
    resp_mat_dict[key] = tmp

# Reshape into single long form dataframe
behav1 = pd.concat(resp_mat_dict, axis=0)
behav1.reset_index(drop=True, inplace=True)

# Add stimulus and subject for indexing
behav1['Tmp'] = list(zip(behav1['CID'], behav1['Experiment'].str[:-1]))
behav1['Stimulus'] = behav1['Tmp'].map(cid_to_stim)
behav1['Subject'] = behav1['Experiment'] + '_' + behav1['ROI#'].astype(str).str.zfill(3)
behav1 = behav1.set_index(['Stimulus', 'Subject']).sort_index()
behav1.drop(['CID', 'Tmp', 'ROI#', 'Experiment'], axis=1, inplace=True)
behav1.head()

# write to disk
molecules.to_csv('molecules.csv')
behav1.to_csv('behavior.csv')
subjects.to_csv('subjects.csv')
stimuli.to_csv('stimuli.csv')
