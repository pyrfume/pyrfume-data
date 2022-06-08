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
raw1 = loadmat('raw/odormappingdata_fullinfo_06_2021.mat')
raw2 = loadmat('raw/odormappingdata_simple_05_2022.mat')
raw3 = loadmat('raw/ID_glom_04_2022.mat')

# raw1 contents:
# 8 structures of form ORdataY_XXXX, 1 for each OB, containing the following:
#  RefImage: 256x256 image of mean baseline fluorescence of the imaged area.
#  RespMatrix: Segmented response matrix for all ROIs (glomeruli) vs. odorants tested. Row is RIO#, column is odor #.
#  OdorList: A list of all of the odors used, in the order that they appear in the RespMatrix. These lists are identical 
#    for all OBs.
#  ROIPos: X (med-lat) and Y (rostr-caud) centroids of each ROI. X from the left OB is negative relative to midline
#  Maps: 256x256 images showing the deltaF response maps for each of the odors (same order as the OdorList and the 
#    RespMatrix).
#  Reg: Info about scaling and registration of OB to reference OB, using midline and caudal sinus as reference landmarks.
#    Pixelsize units are microns/pixel.
#  Metadata: Reference data containing parameters used to generate maps and output values to response matrices.

ref_im, resp_mat, maps = {}, {}, {}
for key in raw.keys():
    if 'ORdata' in key:
        ob = key.split('_')[1]
        tmp = raw[key][0][0]
        ref_im[ob] = np.flip(tmp['RefImage'].T)
        resp_mat[ob] = tmp['RespMatrix']
        maps[ob] = [np.flip(m.T) for m in tmp['Maps'].squeeze()]

# ROIPos from .mat is in table format which can't be loaded by scipy.io.loadmat, must export directly from Matlab to .csv.
# In Matlab use: writetable(ORdataY_XXXX.ROIPos, 'XXXX_ROIPos.csv', 'WriteRowNames', true)
ROIPos = {}
for ob in maps.keys():
    ROIPos[ob] = pd.read_csv('raw/' + ob + '_ROIPos.csv', index_col=0)
        
# Odor lists are identical, only need to grab once. First 185 are the odorants; last 2 are 'empty' (control) and the 
# solvent (triglyceride)
odor_list = [x[0] for x in raw['ORdata1_111L'][0][0]['OdorList'].squeeze()]
odorant_names = [' '.join(x.split(' ')[1:]) for x in odor_list]

# raw2 contents:
# allrespmatrices_dF: Structure array with odorant response matrices and ROI positions for each OB. 
#   Variables are indexed in the following order: (1) omp111L, (2) omp111R, (3) omp112L, (4) omp112R, (5) omp113L, 
#   (6) omp113R, (7) omp114L, (8) omp114R. All odorants are in the same order for all response matrices 
#   (see ‘odornameslist’ variable). ROI positions indicate centroid of each ROI, after visual registration by aligning
#   the midline and caudal sinus. Units are microns, reference (zero) is midline (for ‘X’) and caudal sinus (for ‘Y’). 
#   Note that the X pos for ROIs from the left OB is negative relative to midline.
# allconcs_prep: Calculated final delivered concentration of each odorant, indexed in the same order as response matrices
#   and odornames list, for each mouse, in mols/L. Note that concentrations estimated from vapor pressure, calibrated air
#   dilution, and liquid dilution, assuming ideal behavior.
# odornameslist: List of all 185 odorants (plus two vehicle controls), indexed in same order as odors in response spectrum
#   matrices.

# But only unique data relavtive to "fullinfo" .mat is the odorant concentration data
all_concs = pd.DataFrame(raw2['allconcs_prep'])

# raw3 contents:
# odornameslist: List of all 185 odorants (plus two vehicle controls), indexed in same order as odors in response 
#   spectrum matrices.
# allgloms_odors_sorted_IDglomorder: Structure containing info about identified glomeruli and diagnostic odorants,
#   including additional odorants meeting criteria for sparseness and reliability, but not meeting cirteria for functional
#   identification. Data should match Tables S2 and S3.
# The structure fields are:
#   numgloms50: Number of glomeruli activated by each odorant after thresholding at 50% of max response. Numbers are in 
#     order of preparation (8 OBs total), indexed in the same order as in the 'allrespmatrix' structure from raw2.
#   maxgloms: The identity of the max-activated ROI in each OB, indexed as above.
#   odorname: (self-explanatory; the first two letters are an in-lab identifier and can be ignored).
#   odorindex: Index of odorants in odornames list and in the response spectrum.
#   erroratio: Fraction of ‘errors’ in match between strongest-activated glomerulus in each imaged OB and most-correlated
#     glomerulus response spectrum, of all possible pairwise comparisons.
#   meansparse: <No description provided>
#   meancorr, mediancorr: Mean and median correlation coefficient (Pearson’s r) in response spectrum for the 
#     strongest-activated glomerulus in each imaged OB for a given odorant, calculated across all pairwise comparisons.
#   glomid: ID assigned to each "functionally identified glomeruli".
#   mean_norm_spectrum, median_norm_spectrum: Mean odorant response spectrum for the ID’d glomerulus, normalized within
#     each OB to the response to the given odorant and then averaged across all imaged OBs. Odorants are in the order given
#     in the ‘odornameslist’ variable.
#   mediolateral, rostrocaudal: Positions of max-activated glomeruli in each OB, indexed as in 'maxgloms' field. Units are
#     microns, reference (zero) is midline (for ML units) and caudal sinus (for RC units).
#   MLmean, RCmean: Mean position of each ID’d glomerulus, averaged across each imaged OB after visual registration by 
#     aligning the midline and caudal sinus.

# 80 odorants passed the "highly conservative" criteria for sparseness and reliability, but aach of the above fields has
# length of 85. These last 5 entries consist of 1 blank row  and 4 repeats, and thus will be ignored.
# glomid, mediolaterial rostrocaudal, MLmean, and RCmean only apply to odorants with "functionally identified glomeruli".
# Othersiwe, the odorants here are included as "eliciting consistently sparse activation".

raw3_fields = ['numgloms50', 'maxgloms', 'odorname', 'erroratio', 'meansparse', 'meancorr', 'mediancorr', 'glomid',
               'mean_norm_spectrum', 'median_norm_spectrum', 'mediolateral', 'rostrocaudal', 'MLmean', 'RCmean']

glom_data = {}
for key in raw3_fields:
    if key == 'odorname':
        glom_data[key] = [x[0][0][0] for x in raw3['allgloms_odors_sorted_IDglomorder'][key].squeeze()[:80]]
    elif key in ['numgloms50', 'maxgloms', 'mean_norm_spectrum', 'median_norm_spectrum', 'mediolateral', 'rostrocaudal']:
        glom_data[key] = [x[0] if x.any() else None for x in raw3['allgloms_odors_sorted_IDglomorder'][key].squeeze()[:80]]
    else:
        glom_data[key] = [x[0][0] if x.any() else None for x in raw3['allgloms_odors_sorted_IDglomorder'][key].squeeze()[:80]]

# Display the reference images
fig, axes = plt.subplots(2, 4, figsize=(16, 8))
i = 0
for ob, im in ref_im.items():
    ax = axes.flat[i]
    ax.imshow(im)
    ax.set_title(ob, fontdict={'fontsize':12})
    ax.set_axis_off()
    i += 1

# Display the first few response maps
fig, axes = plt.subplots(8, 4, figsize=(16, 32))
i = 0
for ob, im in maps.items():
    for j in range(0,4):
        ax = axes.flat[i*4 + j]
        ax.imshow(im[j])
        ax.set_title(ob + ': ' + odor_list[j], fontdict={'fontsize':12})
        ax.set_axis_off()
    i += 1

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
           'N,N-dimethyl-2-phenethylamine': 25124,
           '5-methyl heptan-3-one oxime': 90787,
           '2-methoxy-3(5 or 6)-isopropylpyrazine': 71311442,
           'ethyl-2,5-dihydro-4-methylthiazole': 58165301,
           'empty': -1} # this is the control...

for name in man_add:
    cids[name] = man_add[name]

# Get molecule info from cids                   
info_dict = pyrfume.from_cids(list(cids.values()))

# Create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID')
molecules.loc[-1] = [None, None, None, 'empty'] # add the control
molecules.sort_index(inplace=True)
molecules.head()

# Create dataframe from response matrices -> behavior_1.csv
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
behav1 = behav1.set_index(['CID', 'Experiment', 'ROI#']).sort_index()
behav1.head()

# Data for diagnostic odorants and functionally-identifiable glomeruli. Split into two dataframes, one for 
# functionally-identifiable glomeruli (-> behavior_2.csv) and one for odorants that consistently elicit sparse activation
# but did not meet the cirteria for functional identification (-> behavior_3.csv).
data_dict = {}
for key in raw3_fields:
    data_dict[key] = glom_data[key]

# Temporary dataframe to make easier to split
tmp = pd.DataFrame(data_dict)
split_row = tmp[tmp['glomid'].isnull()].index[0]

behav2 = tmp.iloc[:split_row].copy()
behav2['odorname'] = behav2['odorname'].str.split(' ').str[1:].str.join(' ')
behav2['CID'] = behav2['odorname'].map(cids)
behav2.drop(['odorname'], axis=1, inplace=True)
behav2['glomid'] = behav2['glomid'].astype(int)
behav2 = behav2.set_index(['CID']).sort_index()
behav2.head()

behav3 = tmp.iloc[split_row:].copy()
behav3['odorname'] = behav3['odorname'].str.split(' ').str[1:].str.join(' ')
behav3['CID'] = behav3['odorname'].map(cids)
behav3.drop(['odorname', 'glomid', 'mediolateral', 'rostrocaudal', 'MLmean', 'RCmean'], axis=1, inplace=True)
behav3 = behav3.set_index(['CID']).sort_index()
behav3.head()

# Data for subjects.csv
subj_dict = {}
for exp in resp_mat.keys():
    subj_dict[exp] = {'Mouse ID': exp[:-1], 'Hemisphere': exp[-1]}
    
subjects = pd.DataFrame.from_dict(subj_dict, orient='index').sort_index()
subjects.index.names = ['Experiment']
subjects.head()

# Data for stimuli.csv
stimuli = pd.DataFrame(all_concs, copy=True)
stimuli.columns=['111', '112', '113', '114']
stimuli['CID'] = stimuli.index
stimuli['CID'] = stimuli['CID'].map({v: k for v, k in enumerate(odorant_names[:-2])}).map(cids)
stimuli = stimuli.melt(id_vars='CID', var_name='Mouse ID', value_name='Conc. (mols/L)')
stimuli = stimuli.set_index(['CID', 'Mouse ID']).sort_index()
stimuli.head()

# Data for rois.csv
rois = pd.concat(ROIPos, axis=0).reset_index()
rois.rename(columns={'level_0': 'Experiment', 'Row': 'ROI#'}, inplace=True)
rois['ROI#'] = rois['ROI#'].str.split(' ').str[1].astype(int)
rois = rois.set_index(['Experiment', 'ROI#']).sort_index()
rois.head()

# Write images (reference and response maps) to .csv files.
# Becase of total disk space, these will not be on Github, but instead accessible by Dropbox.
im_dir = 'images'
if not os.path.exists(im_dir):
    os.makedirs(im_dir)

mice = list(set([x[:-1] for x in ref_im.keys()]))
for m in mice:
    for i, odor in enumerate(odorant_names):
        file_name = im_dir + '/' + m + '_' + odor.replace(' ', '_') + '.csv'
        pd.DataFrame(maps[m + 'L'][i]).to_csv(file_name, header=None, index=None, float_format="%.8f")

# Dataframe for images.csv, contains file names for maps saved as .csv, indexed by CID
im_data_dict = {}
for exp in resp_mat.keys():
    for odor in odorant_names:
        file_name = im_dir + '/' + exp[:-1] + '_' + odor.replace(' ', '_') + '.csv'
        im_data_dict[(cids[odor], exp)] = file_name
    
images = pd.DataFrame.from_dict(im_data_dict, orient='index', columns=['ImagePath'])
images.index = pd.MultiIndex.from_tuples(images.index)
images.index.names = ['CID', 'Experiment']
images.sort_index(inplace=True)
images.head()

# write to disk
molecules.to_csv('molecules.csv')
images.to_csv('imaging.csv')
behav1.to_csv('behavior_1.csv')
behav2.to_csv('behavior_2.csv')
behav3.to_csv('behavior_3.csv')
subjects.to_csv('subjects.csv')
stimuli.to_csv('stimuli.csv')
rois.to_csv('rois.csv')