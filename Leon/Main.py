#!/usr/bin/env python
# coding: utf-8


import pandas as pd
from pathlib import Path
from PIL import Image
from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import shutil
from collections import OrderedDict
from pyrfume import odorants


# We'll work out of the 'patterns_csv' directory, which has glomerular maps + CAS#s bundled together in individual csv files. First flatten the directory:


saved_path = os.getcwd()
root = '/Users/jcastro/Dropbox/CnG/' # set as needed
patterns_dir = 'Leon/Glomerular archive/patterns_csv' # dir of glomerular maps
path = Path(root + patterns_dir)

# fxn obtained and modified from https://amitd.co/code/python/flatten-a-directory
def flatten(directory):
    for dirpath, _, filenames in os.walk(directory, topdown=False):
        for filename in filenames:
            i = 0
            source = os.path.join(dirpath, filename)
            target = os.path.join(directory, filename)

            while os.path.exists(target):
                i += 1
                file_parts = os.path.splitext(os.path.basename(filename))

                target = os.path.join(
                    directory,
                    file_parts[0] + "_" + str(i) + file_parts[1],
                )

            shutil.move(source, target)

        if dirpath != directory:
            os.rmdir(dirpath)
            
# if dir isn't flattend already:
#flatten(path) 


# Walk through the flattened directory, pulling out image data and the corresponding CAS# from the csv files:



Names = [] # simple file/molecule names (nomenclature idiosyncratic to Leon)
CAS_nums = [] 
Conditions = [] # metadata on experimental conditions (duration, concentration, etc)
Data = [] # list of dfs with the actual data

for file in path.glob('*.csv'):
    df = pd.read_csv(file, index_col=None, header=0)
    # Check to see where data starts
    j=0
    while True:
        if df.iloc[j, 0]=='-100':
            break
        j+=1
    data = df.iloc[j:, :].astype(float)
    data[data<-99] = np.nan
    vals = data.values.ravel()
    assert vals[~np.isnan(vals)].min() > -15
    assert vals[~np.isnan(vals)].max() < 15
    
    name = df.iloc[0, 0]
    cas = df.columns[0]
    condition = df.iloc[1, 0] if j==2 else ''
    
    Names.append(name)
    CAS_nums.append(cas)
    Conditions.append(condition) 
    Data.append(data)


# Clean up the molecule lists, and write csvs of glomerular maps to a new tmp directory. 
# For duplicate maps, just name them as X_replicate_1, X_replicate_2, etc.  



# first, remove list indices where the CAS# is = 0 (these probably correspond to mixtures, or things like 'coffee', etc)
zeros_idx_list = sorted([ i for i in range(len(CAS_nums)) if CAS_nums[i] == '0' ], reverse=True)

for idx in zeros_idx_list:
    Names.pop(idx)
    CAS_nums.pop(idx)
    Conditions.pop(idx)
    Data.pop(idx)
    

# extract unique molecule names, and find replicates     
tmp_path = '/Users/jcastro/Dropbox/CnG/Leon/CSVs_for_Pyrfume/'
isExist = os.path.exists(tmp_path)

if not isExist:
    os.mkdir(tmp_path)

cas_unique = np.unique(CAS_nums)
j=0
behav = []
rep_string = []



for m in range(len(cas_unique)):
    indices = [i for i, x in enumerate(CAS_nums) if x == cas_unique[m]]
    # only one map corresponding to the odorant
    if len(indices) == 1:
        molName = Names[indices[0]] + '.csv'
        Data[indices[0]].to_csv(tmp_path + molName, index=False,header=False)
        behav.append(molName)
        
    # multiple maps corresponding to the odorant    
    else:
        j = 0
        rep_string = ''
        for id in indices:
            j = j+1
            molName = Names[indices[0]] + str('_replicate_'+ str(j) + '.csv')
            rep_string += (molName + ';')
            Data[id].to_csv(tmp_path + molName, index=False,header=False)
        behav.append(rep_string) 


# For whatever reason, one CAS# (for Methyl 3-aminobenzoate) is just formatted really strangely. Clean up this one pathological case: 



# fix a bad CAS#:
cas_array = cas_unique.tolist()
idx = cas_array.index('10/9/4518')

cas_unique[idx] = '4518-10-9'


# Fetch CIDs from CAS#s, perform standardization for Pyrfume: 


molecules = OrderedDict()
molecules = odorants.get_cids(cas_unique)


# grab other identifiers for each molecule (SMILES, IUPAC, etc):
mols = odorants.from_cids(list(molecules.values()))
molecules_df = pd.DataFrame(mols)
molecules_df.set_index('CID', inplace = True)


# Write the 'molecule' and 'behavior' files for pyrfume



tmp_path = '/Users/jcastro/Dropbox/CnG/Leon/CSVs_for_Pyrfume/'

# Molecules
molecules_df.to_csv('molecules_fin.csv')

# Behavior
dd = pd.DataFrame(behav)
dd.to_csv(tmp_path + 'behavior_fin.csv',index=False)

