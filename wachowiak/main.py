#!/usr/bin/env python
# coding: utf-8

get_ipython().run_line_magic('matplotlib', 'inline')
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyrfume.odorants import get_cids, from_cids
from scipy.io import loadmat
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style('dark')

# load raw data from .mat file and convert to dict
raw = loadmat('omp114_Lbulb_allglomsORdata_final_registered_checked0423.mat')
raw = raw['ORdata'][0][0]

# raw contents:
# RefImage: Just an image of the resting fluorescence of the prep.
# RespMatrix: A matrix of the response values of every ROI (glomerulus) for every odorant tested.  Row is RIO#, column is odor #.
# OdorList: A list of all of the odors used, in the order that they appear in the RespMatrix.
# ROIPos: A table of the X and Y coordinates for the centroids of each of the ROIs.
# Maps: An array of 256x256 pixel images showing the delta F response maps for each of the odors (same order as the OdorList and the RespMatrix).
# Reg: You can ignore this for now - used to registering maps from different preparations.
# Metadata: Also, can ignore - just has internal ref info for us regarding the experiment.

d = {}
keys = ['RefImage', 'RespMatrix', 'OdorList', 'Maps', 'Reg', 'MetaData']
for key in keys:
    d[key] = raw[key].squeeze()

d['OdorList'] = [x[0] for x in d['OdorList']]

# ROIPos from .mat is in table format which can't be loaded by scipy.io.loadmat
# From Matlab, ROIPos was exported to .csv so that it can be opened by Pandas

ROIPos = pd.read_csv('ROIPos.csv', index_col=0)

# display the reference image
plt.imshow(np.flip(d['RefImage'].T));
plt.axis('off');

# display the first few response maps
fig, axes = plt.subplots(5, 5, figsize=(12, 12))
for i, ax in enumerate(axes.flat):
    ax.imshow(np.flip(d['Maps'][i].T))
    ax.set_title(d['OdorList'][i], fontdict={'fontsize':12})
    ax.set_axis_off()

# get cids
odorant_names = [' '.join(x.split(' ')[1:]) for x in d['OdorList']]
cids = get_cids(odorant_names)

# manually add the missing CIDs
man_add = {'(+-)-geosmin': 15559490,
           '5-methyl heptan-3-one oxime': 90787,
           '2-methoxy-3(5 or 6)-isopropylpyrazine': 71311442,
           'ethyl-2,5-dihydro-4-methylthiazole': 58165301,
           'empty': -1} # presumably this is the control...

for name in man_add:
    cids[name] = man_add[name]

# get molecule info from cids                   
info_dict = from_cids(list(cids.values()))

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# Create dataframe from RespMatrix
resp_mat = pd.DataFrame(d['RespMatrix'], index=ROIPos.index, columns=d['OdorList'])
resp_mat.head()

# Heatmap of response matrix
sns.set_style("ticks", {"xtick.bottom": True, "ytick.left": True})
plt.figure(figsize=(30, 5))
ax = sns.heatmap(resp_mat, xticklabels=resp_mat.columns, yticklabels=20, cbar=None, vmax=50)
ax.set_xticklabels(ax.get_xticklabels(), fontsize=7, ha='left');
ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, va='center');
plt.xlabel('Odorant')
plt.ylabel('Glomerulus (ROI)');

# write response maps to separate .csv files
imDir = 'images'
if not os.path.exists(imDir):
    os.makedirs(imDir)
    
for i, im in enumerate(d['Maps']):
    fileName = imDir + '\\' + str(d['OdorList'][i]).replace(' ', '_') + '.csv'
    pd.DataFrame(np.flip(im.T)).to_csv(fileName, header=None, index=None, float_format="%.8f")

# dataframe for images.csv, contains file names for maps saved as .csv, indexed by CID
data_dict = {}

for i, odor in enumerate(odorant_names):
    fileName = imDir + '\\' + str(d['OdorList'][i]).replace(' ', '_') + '.csv'
    data_dict[cids[odor]] = fileName
    
images = pd.DataFrame.from_dict(data_dict, orient='index', columns=['ImagePath']).sort_index()
images.index.name = 'CID'
images.head()

# write to disk
molecules.to_csv('molecules.csv')
images.to_csv('images.csv')
resp_mat.to_csv('resp_mat.csv')