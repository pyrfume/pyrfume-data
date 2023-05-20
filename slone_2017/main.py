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

import pandas as pd
import numpy as np
import pyrfume

# +
# Dataset S1: Average SSR and EAG data
excel_settings = [ # (sheet name, header row, adjustment type)
    ('SSR (Pentane Adjusted)', 3, 'SSR_Pentane'),
    ('SSR (Adjusted to "No UAS")', 2, 'SSR_No-UAS'),
    ('EAG (Paraffin Oil Adjusted)', 2, 'EAG_Paraffin'),
    ('EAG (Adjusted to "No UAS")', 2, 'EAG_No-UAS')    
]

ds1 = {}
for sheet, header, adjustment in excel_settings:
    df = pd.read_excel('pnas.1704647114.sd01.xlsx', sheet_name=sheet, header=header, index_col=0, skipfooter=1)
    df.index.name = 'Stimulus'
    df = df.iloc[:, :-1]
    df = pd.melt(df, var_name='Subject', value_name=adjustment, ignore_index=False)
    df.Subject = df.Subject.str.replace(' ', '')
    ds1[adjustment] = df.reset_index().set_index(['Stimulus', 'Subject']).sort_index()

# +
# SSR data -> behavior_1.csv
ssr = ds1['SSR_Pentane'].join(ds1['SSR_No-UAS'])
ssr = pd.melt(ssr, var_name='Normalization', value_name='Ave SSR (DeltaSpikes/s)', ignore_index=False)
ssr['Normalization'] = ssr['Normalization'].str.split('_').str[-1]

ssr.head()

# +
# EAG data -> behavior_2.csv
eag = ds1['EAG_Paraffin'].join(ds1['EAG_No-UAS'])
eag = pd.melt(eag, var_name='Normalization', value_name='Ave EAG Response', ignore_index=False)
eag['Normalization'] = eag['Normalization'].str.split('_').str[-1]

eag.head()

# +
# Dataset S2: PCA components -> behavior_3.csv
ds2 = pd.read_excel('pnas.1704647114.sd02.xlsx', sheet_name=None, header=0, index_col=0)

pca = pd.concat([ds2['SSR'], ds2['EAG']], axis=0)
pca.index.name = 'Stimulus'

pca.head()

# +
# Load odorants table
odorants = pd.read_excel('Table_S2.xlsx')

odorants.head()
# -

# Get CIDs from names
cids_from_name = pyrfume.get_cids(odorants['Full name'].tolist(), kind='name')

# +
# Stimuli
stimuli = odorants.copy().set_index('Abbreviation').sort_index()
stimuli.index.name = 'Stimulus'
stimuli['CID'] = stimuli['Full name'].map(cids_from_name)
# No CID for 'Henkel 100' or 'C7 to C40 Mixture'
stimuli.replace({0: None, 'Not available': None}, inplace=True)

stimuli.head()
# -

# Molecule info
cids = stimuli.CID.dropna().astype(int).to_list()
molecules = pd.DataFrame(pyrfume.from_cids(cids)).set_index('CID').sort_index()
molecules.head()

# +
# Subject info is the same for all expriments; get from 1st sheet from Dataset 1
subjects = pd.read_excel(
    'pnas.1704647114.sd01.xlsx',
    sheet_name='SSR (Pentane Adjusted)',
    header=None,
    index_col=0,
    nrows=3
).T

# Need to manually enter 'Enrichment' since values are color coded (workers = orange dots; males = green dots)
subjects['Enrichment'] = [None] * 2 + ['worker'] * 5 + [None] * 2 + ['male'] + [None] * 4 + ['male'] * 3 + [None] + \
    ['male'] * 3 + ['worker'] + [None] * 2 + ['male', None]

subjects['Subject'] = subjects['OR Line'].str.replace(' ', '')
subjects = subjects.set_index('Subject')
subjects.loc['NoUAS', 'Subfamily'] = 'Control'

subjects.head()
# -

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
ssr.to_csv('behavior_1.csv')
eag.to_csv('behavior_2.csv')
pca.to_csv('behavior_3.csv')
subjects.to_csv('subjects.csv')
