#!/usr/bin/env python
# coding: utf-8


import os
import scipy.io as sio
import pandas as pd
from pyrfume import odorants
import numpy as np
from collections import OrderedDict


# ### Description of Data: 
# 
# There are 4 datasets here: 
# 
# 1. **Glomerular data:** glomerular deltaF/F responses organized into 57 (odorants) x N (glomeruli) matrices. Data are from 5 animals, for both left and right hemibulbs (i.e. 10 x 57 x N). N (the # of glomeruli) is variable for each (hemibulb, animal) pair. 
# 
# 
# 2. **'Mitral 55' data:** cellular deltaF/F responses organized into 55 (odorants) x N (cells) matrices. There are data for both 'high' and 'low' concentrations of odorant, for 7 and 6 fields of view, respectively. Similar to above, N is the number of cells, and varies across concentration and FOV.
# 
# 
# 3. **'Mitral 33' data:** cellular deltaF/F responses organized into 33 (odorants) x N (cells) matrices. Only one concentration, but 6 fields of view (6 total conditions). N = # cells, and varies across the 6 conditions. 
# 
# 
# 4. **Tufted data:** cellular deltaF/F rsponses, organized into 55 (odorants) x N (cells) matrices. 2 concentrations and 3 fields of view (6 total conditions). N = # cells, and varies across the 6 conditions. 
# 
# 

# ### Glomerular data: Delta F/F responses

# First grab the CIDs for the glomerular data (55 odorants)


saved_path = os.getcwd()
root = '/Users/jcastro/Dropbox/CnG/Chae et al 2019/' # set as needed
glom_data_path = 'raw/glomerular_data/glomerular_responses/glomerular_responses__main'
mol_names_path = root + 'raw/supplementary table 1 odors with cas and cid.xlsx'

os.chdir(root + glom_data_path)

glom_CIDs = pd.read_excel(mol_names_path, header = 0, sheet_name='Glomerular odor panel',
                              usecols = 'C', squeeze=True)

# putting  a dummy value in for the nan, so that the function odorants.from_cids doesn't throw an error
glom_CIDs = glom_CIDs.replace(np.nan, 123456).values.tolist()
glom_CIDs = list(map(int,glom_CIDs))

print(root + glom_data_path)
os.chdir(saved_path)


# Grab all the dF/F data for the glomerular experiments, and assemble into long format


# walk through the directories and grab the glomerular response data
i=0

for root, dirs, files in os.walk(root + glom_data_path):
    for filename in files:
        res = root.split('/')
        hemibulb, animal = res[-1], res[-2]
        if filename.endswith('responses.mat'):
            thedata = sio.loadmat(root + '/'+ filename)
            hemibulb_string = hemibulb + '_ROI_responses'
            thedata = thedata[hemibulb_string]
            df = pd.DataFrame(thedata)
            df.insert(0, 'animal', animal[-1]) # column w/ experimental subject 
            df.insert(0, 'hemibulb', hemibulb.split('_')[0]) # column indicating which hemisphere
            df.insert(0,'CIDs',glom_CIDs)
            # melt into long format:
            _, cols = thedata.shape
            df = pd.melt(df, id_vars=['CIDs', 'hemibulb', 'animal'], value_vars=list(range(cols)))
            df.rename(columns={'variable': 'glom', 'value' : 'delF'}, inplace=True)
            df.sort_values(by=['CIDs'], inplace=True)
            
            if i == 0: 
                df_cat_glom = df
            else:
                df_cat_glom = pd.concat([df_cat_glom, df])
            i = i + 1

df_cat_glom.sort_values(by=['CIDs', 'animal', 'glom'], inplace=True)
df_cat_glom.shape


# Quick and dirty function to extract data by animal and hemibulb (the orig. format)



def extractGlomData(df,LorR,animal_num):
    cond1 = df.hemibulb == LorR
    cond2 = df.animal == animal_num
    all_cond = cond1 & cond2
    df_out = df[all_cond].pivot(index = 'CIDs', columns=['glom'], values = 'delF')
    return df_out

# example: 
df_out = extractGlomData(df_cat_glom,'right','1') # return the df for right hemibulb of animal #1
print(df_out.shape) # for ('right', animal # 1) should return 57, 116
df_out.head()


# ### Fetch the 'Mitral 55', 'Mitral 33', and 'tufted' data: Delta F/F responses

# The 3 datasets are all sitting in the same directory, so we'll grab each of them in the same loop and 
# save them as dicts. 



saved_path = os.getcwd()
root = '/Users/jcastro/Dropbox/CnG/Chae et al 2019/' # set as needed
mitral = 'raw/MT_cell/'

#mitral_cells/'
os.chdir(root + mitral)

mitral_33 = OrderedDict()
mitral_55 = OrderedDict()
tufted = OrderedDict()

# walk through the directories and grab the glomerular response data
for root, dirs, files in os.walk(root + mitral):
    for filename in files:
        res = root.split('/')
        experiment_type = res[-1] # tufted, or mitral1 (33 odorants), or mitral2 (55 odorants)
        
        if filename.endswith('.mat'):
            namestring = filename.split('.')[-2]
            thedata = sio.loadmat(root + '/'+ filename) 
            
            if experiment_type == 'tufted_cells':
                tufted[namestring] = thedata['R']
            if experiment_type == '33_odors':
                mitral_33[namestring] = thedata['R']
            if experiment_type == '55_odors':
                mitral_55[namestring] = thedata['R']

mitral33_mol_list = pd.read_excel('MT_cells.xlsx', sheet_name='odor list_33 odors', header=0, usecols = 'B')
mitral55_mol_list = pd.read_excel('MT_cells.xlsx', sheet_name='odor list_55 odors', header=0, usecols = 'B')
tufted_mol_list = '' # fix this... 

os.chdir(saved_path)


# ### 'Mitral 55', 'Mitral 33'molecule lists

# Generate the list of Mitral55 odorants:



mitral55_mol_list
molecules_55 = odorants.get_cids(list(mitral55_mol_list['odor name']))



# Fix the pathological cases where pyrfume didn't find a CID: 

molecules_55['ethyl 3-mercapto propionate'] = 87475436  
molecules_55['fenchone (-)'] = 14525  
molecules_55['citral cis+trans'] = 638011  
molecules_55['3-acetal 2,5 dimethyl furan'] = 61527
molecules_55['phenyl ethyl acetate'] = 854188  
molecules_55['carvyl acetate (-)'] = 11964245  
molecules_55['trans 2-hexenol'] = 5318042

mitral55_cids = list(molecules_55.values())


# Generate the list of Mitral33 odorants: 


mitral33_mol_list
molecules_33 = odorants.get_cids(list(mitral33_mol_list['odor name']))



# Fix the pathological cases where pyrfume didn't find a CID: 

molecules_33['ethyl 3-mercapto propionate'] = 87475436  
# going on the assumption that this is mislabeled, and should be 'isopropyl', not 'siopropyl'
molecules_33['4-siopropyl benzaldehyde'] = 326

mitral33_cids = list(molecules_33.values())


# ### Tidy the 'mitral 55', 'mitral 33', and 'tufted' data: 

# Mitral 55 data


keys = mitral_55.keys()
i = 0

for key in keys:
    HorL = key[1] #l (low) or h (high) concentration
    fov = key[-1] # field of view
    thedata = mitral_55[key]
    df = pd.DataFrame(thedata)
    df.insert(0, 'High or Low Conc.', HorL) # column with conc. level (h (high) or l (low))
    df.insert(0, 'FOV', fov) # column with FOV number (see orig paper)
    df.insert(0,'CIDs',mitral55_cids)
    
    # melt into long format:
    _, cols = thedata.shape
    df = pd.melt(df, id_vars=['CIDs', 'FOV', 'High or Low Conc.'], value_vars=list(range(cols)))
    df.rename(columns={'variable': 'cell', 'value' : 'delF'}, inplace=True)
    df.sort_values(by=['CIDs'], inplace=True)
    
    if i == 0: 
        df_cat_mitral55 = df
    else:
        df_cat_mitral55 = pd.concat([df_cat_mitral55, df])
    i = i + 1

df_cat_mitral55.sort_values(by=['CIDs', 'FOV', 'cell'], inplace=True)


# Mitral 33 data



keys = mitral_33.keys()
i = 0

for key in keys:
    fov = key[-1] # field of view
    thedata = mitral_33[key]
    df = pd.DataFrame(thedata)
    df.insert(0, 'FOV', fov) # column with FOV number (see orig paper)
    df.insert(0,'CIDs',mitral33_cids)
    
    # melt into long format:
    _, cols = thedata.shape
    df = pd.melt(df, id_vars=['CIDs', 'FOV'], value_vars=list(range(cols)))
    df.rename(columns={'variable': 'cell', 'value' : 'delF'}, inplace=True)
    df.sort_values(by=['CIDs'], inplace=True)
    
    if i == 0: 
        df_cat_mitral33 = df
    else:
        df_cat_mitral33 = pd.concat([df_cat_mitral33, df])
    i = i + 1

df_cat_mitral33.sort_values(by=['CIDs', 'FOV', 'cell'], inplace=True)


# Tufted data




#Tufted data: cellular deltaF/F rsponses, organized into 55 (odorants) x N (cells) matrices. 
#2 concentrations and 3 fields of view (6 total conditions). N = # cells, and varies across the 6 conditions.

keys = tufted.keys()
i = 0

for key in keys:
    HorL = key[1] #l (low) or h (high) concentration
    fov = key[-1] # field of view
    thedata = tufted[key]
    df = pd.DataFrame(thedata)
    df.insert(0, 'High or Low Conc.', HorL) # column with conc. level (h (high) or l (low))
    df.insert(0, 'FOV', fov) # column with FOV number (see orig paper)
    df.insert(0,'CIDs',mitral55_cids)
    
    # melt into long format:
    _, cols = thedata.shape
    df = pd.melt(df, id_vars=['CIDs', 'High or Low Conc.', 'FOV'], value_vars=list(range(cols)))
    df.rename(columns={'variable': 'cell', 'value' : 'delF'}, inplace=True)
    df.sort_values(by=['CIDs', 'cell'], inplace=True)
    
    if i == 0: 
        df_cat_tufted = df
    else:
        df_cat_tufted = pd.concat([df_cat_tufted, df])
    i = i + 1

df_cat_tufted.sort_values(by=['CIDs', 'FOV', 'High or Low Conc.', 'cell'], inplace=True)


# ### Merging Molecule lists




# generate merged molecule lists:
glom_cids = glom_CIDs
mitral55_cids = list(molecules_55.values())
mitral33_cids = list(molecules_33.values())
all_mols = list(set(glom_cids + mitral55_cids + mitral33_cids))





# fetch IUPAC, SMILES, etc. info
mols = odorants.from_cids(all_mols)





mols_df = pd.DataFrame(mols) # Still need to fix 'coffee'
mols_df.head(6)


# Recall that CID='123456' was a dummy value for 'coffee', which
# doesn't actually have a CID. So we need to adjust the row w/ that
# CID entry (row # 36): 




mols_df.iloc[36] = [-1, '', '', '', 'coffee']
# write to disk
mols_df.sort_values(by='CID', inplace=True)


# While we're attending to this, let's also fix the CID entries from the glomerular data that corresponde to coffee, and have the incorrect dummy value of 123456. The correct CID should be -1:




coffee_idx = (df_cat_glom['CIDs'] == 123456)
df_cat_glom.CIDs[coffee_idx] = -1
df_cat_glom.sort_values(by=['CIDs', 'animal', 'glom'], inplace=True)


# ### Write behavior and molecule files to disk:



df_cat_mitral55.to_csv('behavior_3.csv', index=False)
df_cat_mitral33.to_csv('behavior_2.csv', index=False)
df_cat_tufted.to_csv('behavior_4.csv', index=False)
df_cat_glom.to_csv('behavior_1.csv', index=False)

mols_df.to_csv('molecules.csv',index=False)

