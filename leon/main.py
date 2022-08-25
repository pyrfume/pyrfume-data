#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import glob
import pyrfume

# Load the data from the glomerular archive, and put into a DF
data = {}
for file in glob.glob('raw/*.csv'):
    df = pd.read_csv(file, index_col=None, header=None, skip_blank_lines=False)
   
    # Parse file contents
    cas = df.iloc[0,0] # CAS in 1st columns, 1st row
    name = df.iloc[1, 0] # molecule name in 1st column, 2nd row
    conditions = df.iloc[2,0] # experimental conditions in 1st column, 3rd row
    
    df = df.iloc[3:, :].astype(float) # Data starts on 4th row (index from 0)
    df[df == -100] = np.nan
    
    data[os.path.split(file)[1]] = [name, cas, conditions, df]

df = pd.DataFrame.from_dict(data, orient='index', columns=['Name', 'CAS', 'Conditions', 'Data'])
df.index.name = 'SourceFile'
df.reset_index(inplace=True)

# Fix a bad CAS and replace '0' with NaN
df.loc[df.CAS == '10/9/4518', 'CAS'] = '4518-10-9'
df.CAS.replace('0', np.nan, inplace=True)

df.head()

# Fetch CIDs from CAS and add to the dataframe
cids = pyrfume.get_cids(df.CAS.dropna().to_list())

# Add CIDs to dataframe and assign values for those missing
df['CID'] = df.CAS.map(cids)
df.loc[df.CID.isna(), 'CID'] = 0

# Use negative numbers for substances
names_missing_cids = df[(df.CID == 0) & (~df.Name.duplicated())].drop('CID', axis=1)
names_missing_cids['CID'] = [-i for i in range(1, len(names_missing_cids)+1)]

df.loc[df.CID == 0, 'CID'] = df.loc[df.CID == 0, 'Name']     .map(dict(zip(names_missing_cids.Name, names_missing_cids.CID)))
df.CID = df.CID.astype(int)

# Add column to track replicate #
df['Rep'] = df.groupby('CID').cumcount()

# Create stimulus ID from CID and Rep#
df['Stimulus'] = df.CID.astype(str) + '_' + df.Rep.astype(str)

# Create output file name for saving activity maps to .csv
df['Outfile'] = 'csvs/' + df.Stimulus + '.csv'

df = df.set_index('Stimulus').sort_index()
df.head()

# Write csv files (glomerular activity maps) to disk if they don't already exist
if not os.path.isdir('csvs'):
    os.makedirs('csvs')
    
    for row in df[['Outfile', 'Data']].itertuples():
        row.Data.to_csv(row.Outfile, header=False, index=False)
else:
    print('CSV files already exist.')

# Get info for molecules.csv
valid_cids = df.loc[df.CID > 0, 'CID'].to_list()
molecules = pd.DataFrame(pyrfume.from_cids(valid_cids)).set_index('CID').sort_index()

molecules.head()

# Create dataframe for stimuli.csv
stimuli = df[['CID', 'Rep', 'Name', 'Conditions', 'SourceFile']].copy()
stimuli.head()


# Dataframe for imaging.csv
imaging = df['Outfile'].copy().to_frame()
imaging.columns = ['Activity Map Path']
imaging.head()

# Write to disk
molecules.to_csv('molecules.csv')
imaging.to_csv('behavior_1.csv')
stimuli.to_csv('stimuli.csv')

