#!/usr/bin/env python
# coding: utf-8

# Description of Chae 2019 ata: 
#
# 1. Glomerular data: glomerular deltaF/F responses organized into 57 (odorants) x N (glomeruli) matrices. Data are from 5
# animals, for both left and right hemibulbs (i.e. 10 x 57 x N). N (the # of glomeruli) is variable for each (hemibulb, animal)
# pair. 
#
# 2. 'Mitral 55' data: cellular deltaF/F responses organized into 55 (odorants) x N (cells) matrices. There are data for both
# 'high' and 'low' concentrations of odorant, for 7 and 6 fields of view, respectively. Similar to above, N is the number of
# cells, and varies across concentration and FOV.
#
# 3. 'Mitral 33' data: cellular deltaF/F responses organized into 33 (odorants) x N (cells) matrices. Only one concentration, 
# but 6 fields of view (6 total conditions). N = # cells, and varies across the 6 conditions. 
#
# 4. Tufted data: cellular deltaF/F rsponses, organized into 55 (odorants) x N (cells) matrices. 2 concentrations and 3 fields
# of view (6 total conditions). N = # cells, and varies across the 6 conditions. 

import os
import scipy.io as sio
import pandas as pd
import pyrfume
import numpy as np
from collections import OrderedDict

# Data directories
raw_root = 'raw'
mol_names_path = os.path.join(raw_root, 'supplementary_table_1_odors_with_cas_and_cid.xlsx')
glom_data_path = os.path.join(raw_root, 'glomerular_data', 'glomerular_responses', 'glomerular_responses__main')
mitral_data_path = os.path.join(raw_root, 'MT_cell')

# CIDs for each experiment
glom_cids = pd.read_excel(mol_names_path, header = 0, sheet_name='Glomerular odor panel')['CID']
glom_cids = glom_cids.replace(np.nan, -1).astype(int).values.tolist() # Use -1 for coffee's missing CID

mt_cell_file = os.path.join(mitral_data_path, 'MT_cells.xlsx')
mol_list_33 = pd.read_excel(mt_cell_file, sheet_name='odor list_33 odors', header=0)['odor name'].to_list()
mol_list_55 = pd.read_excel(mt_cell_file, sheet_name='odor list_55 odors', header=0)['odor name'].to_list() 

molecules_55_name_to_cid = odorants.get_cids(mol_list_55, kind='name')

# Fix the pathological cases where pyrfume didn't find a CID: 
molecules_55_name_to_cid['ethyl 3-mercapto propionate'] = 87475436  
molecules_55_name_to_cid['fenchone (-)'] = 14525  
molecules_55_name_to_cid['citral cis+trans'] = 638011  
molecules_55_name_to_cid['3-acetal 2,5 dimethyl furan'] = 61527
molecules_55_name_to_cid['phenyl ethyl acetate'] = 854188  
molecules_55_name_to_cid['carvyl acetate (-)'] = 11964245  
molecules_55_name_to_cid['trans 2-hexenol'] = 5318042
molecules_55_name_to_cid['coffee'] = -1

cids_55 = list(molecules_55_name_to_cid.values())

molecules_33_name_to_cid = odorants.get_cids(mol_list_33, kind='name')

# Fix the pathological cases where pyrfume didn't find a CID: 
molecules_33_name_to_cid['ethyl 3-mercapto propionate'] = 87475436 

# Going on the assumption that this is mislabeled, and should be 'isopropyl', not 'siopropyl'
molecules_33_name_to_cid['4-siopropyl benzaldehyde'] = 326

cids_33 = list(molecules_33_name_to_cid.values())

# Merge CID lists and get info for molecules.csv
all_cids = list(set(glom_cids + cids_55 +cids_33))

molecules = pd.DataFrame(pyrfume.from_cids(all_cids)).set_index('CID').sort_index()

molecules.head()

# Grab all the dF/F data for the glomerular experiments, and assemble into long format
# Walk through the directories and grab the glomerular response data

glom_dict = {}
for root, dirs, files in os.walk(glom_data_path):
    for filename in files:
        animal, hemibulb = os.path.normpath(root).split(os.sep)[-2:]
        if filename.endswith('responses.mat'):
            the_data = sio.loadmat(os.path.join(root, filename))
            df = pd.DataFrame(the_data[hemibulb + '_ROI_responses'])
            df.insert(0, 'Animal', int(animal[-1])) # column w/ experimental subject 
            df.insert(0, 'Hemibulb', hemibulb.split('_')[0]) # column indicating which hemisphere
            df.insert(0, 'CID', glom_CIDs)
            
            # melt into long format:
            df = pd.melt(df, id_vars=['CID', 'Hemibulb', 'Animal'], var_name='Glom', value_name='DeltaF/F')
            
            glom_dict[animal + '_' + hemibulb] = df

glom_data = pd.concat(glom_dict, axis=0)

glom_data['Stimulus'] = 'G_' + glom_data.CID.astype(str)
glom_data['Subject'] = 'G_' + glom_data[['Animal', 'Hemibulb']].astype(str).agg('_'.join, axis=1) + '_'     + glom_data.Glom.astype(str).str.zfill(2)

glom_data = glom_data.reset_index(drop=True).set_index(['Stimulus', 'Subject']).sort_index()
glom_data.head()

# Function to extract data by animal and hemibulb (the orig. format)
def extractGlomData(df, LorR, animal_num):
    return df[(df.Hemibulb == LorR) & (df.Animal == animal_num)].pivot(index = 'CID', columns=['Glom'], values = 'DeltaF/F')

# example: 
df_out = extractGlomData(glom_data, 'right', 1) # return the df for right hemibulb of animal #1
print(df_out.shape) # for ('right', animal # 1) should return 57, 116
# df_out.head()

# Fetch the 'Mitral 55', 'Mitral 33', and 'tufted' data: Delta F/F responses
# The 3 datasets are all same parent directory, so we'll grab each of them in the same loop and save them as dicts
mitral_33 = OrderedDict()
mitral_55 = OrderedDict()
tufted = OrderedDict()

# Walk through the directories and grab the glomerular response data
for root, dirs, files in os.walk(mitral_data_path):
    for filename in files:
        if filename.endswith('.mat'):
            the_data = sio.loadmat(os.path.join(root, filename)) 
            
            experiment_type = os.path.normpath(root).split(os.sep)[-1] # tufted, mitral1 (33 odorants), or mitral2 (55 odorants)  
            namestring = filename.split('.')[0]
    
            if experiment_type == 'tufted_cells':
                tufted[namestring] = the_data['R']
            if experiment_type == '33_odors':
                mitral_33[namestring] = the_data['R']
            if experiment_type == '55_odors':
                mitral_55[namestring] = the_data['R']

# Mitral 55 data
mitral55_dict = {}
for key in mitral_55.keys():
    df = pd.DataFrame(mitral_55[key])
    df.insert(0, 'High/Low Conc.',key[1]) # column with conc. level (h (high) or l (low))
    df.insert(0, 'FOV', key[-1]) # column with FOV number (see orig paper)
    df.insert(0, 'CID', cids_55)
    
    # melt into long format:
    df = pd.melt(df, id_vars=['CID', 'FOV', 'High/Low Conc.'], var_name='Cell', value_name='DeltaF/F')
    
    mitral55_dict[key] = df
    
mitral55_data = pd.concat(mitral55_dict, axis=0)

mitral55_data['Stimulus'] = 'M55_' + mitral55_data.CID.astype(str) + '_' + mitral55_data['High/Low Conc.']
mitral55_data['Subject'] = 'M55_' + mitral55_data[['Cell', 'FOV']].astype(str).agg('_'.join, axis=1)

mitral55_data = mitral55_data.reset_index(drop=True).set_index(['Stimulus', 'Subject']).sort_index()
mitral55_data.head()

# Mitral 33 data
mitral33_dict = {}
for key in mitral_33.keys():
    df = pd.DataFrame(mitral_33[key])
    df.insert(0, 'FOV', key[-1]) # column with FOV number (see orig paper)
    df.insert(0,'CID',  cids_33)
    
    # melt into long format:
    df = pd.melt(df, id_vars=['CID', 'FOV'], var_name='Cell', value_name='DeltaF/F')
        
    mitral33_dict[key] = df

mitral33_data = pd.concat(mitral33_dict, axis=0)

mitral33_data['Stimulus'] = 'M33_' + mitral33_data.CID.astype(str)
mitral33_data['Subject'] = 'M33_' + mitral33_data[['Cell', 'FOV']].astype(str).agg('_'.join, axis=1)

mitral33_data = mitral33_data.reset_index(drop=True).set_index(['Stimulus', 'Subject']).sort_index()
mitral33_data.head()

# Tufted data
tufted_dict = {}
for key in tufted.keys():
    df = pd.DataFrame(tufted[key])
    df.insert(0, 'High/Low Conc.',key[1]) # column with conc. level (h (high) or l (low))
    df.insert(0, 'FOV', key[-1]) # column with FOV number (see orig paper)
    df.insert(0, 'CID', cids_55)
    
    # melt into long format:
    df = pd.melt(df, id_vars=['CID', 'FOV', 'High/Low Conc.'], var_name='Cell', value_name='DeltaF/F')
    
    tufted_dict[key] = df
    
tufted_data = pd.concat(tufted_dict, axis=0)

tufted_data['Stimulus'] = 'T_' + tufted_data.CID.astype(str) + '_' + tufted_data['High/Low Conc.']
tufted_data['Subject'] = 'T_' + tufted_data[['Cell', 'FOV']].astype(str).agg('_'.join, axis=1)

tufted_data = tufted_data.reset_index(drop=True).set_index(['Stimulus', 'Subject']).sort_index()
tufted_data.head()

# Create dataframe for stimuli.csv
stimuli = pd.concat([glom_data.reset_index()[['Stimulus', 'CID']].copy(),
                    mitral33_data.reset_index()[['Stimulus', 'CID']].copy(),
                    mitral55_data.reset_index()[['Stimulus', 'CID', 'High/Low Conc.']].copy(),
                    tufted_data.reset_index()[['Stimulus', 'CID', 'High/Low Conc.']].copy()], axis=0)

stimuli = stimuli.set_index('Stimulus').sort_index()
stimuli = stimuli[~stimuli.index.duplicated()]
stimuli.head()

# Create dataframe for subjects.csv
subjects = pd.concat([glom_data.reset_index()[['Subject', 'Animal', 'Hemibulb', 'Glom']].copy(),
                     mitral33_data.reset_index()[['Subject', 'Cell', 'FOV']].copy(),
                     mitral55_data.reset_index()[['Subject', 'Cell', 'FOV']].copy(),
                     tufted_data.reset_index()[['Subject', 'Cell', 'FOV']].copy()], axis=0)

subjects = subjects.set_index('Subject').sort_index()
subjects = subjects[~subjects.index.duplicated()]
subjects.head()

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
subjects.to_csv('subjects.csv')
glom_data.drop(['CID', 'Hemibulb', 'Animal', 'Glom'], axis=1).to_csv('behavior_1.csv')
mitral33_data.drop(['CID', 'FOV', 'Cell'], axis=1).to_csv('behavior_2.csv')
mitral55_data.drop(['CID', 'FOV', 'High/Low Conc.', 'Cell'], axis=1).to_csv('behavior_3.csv')
tufted_data.drop(['CID', 'FOV', 'High/Low Conc.', 'Cell'], axis=1).to_csv('behavior_4.csv')

