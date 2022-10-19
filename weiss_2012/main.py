#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import tabula
import pyrfume

# Load odorants used; read in Tables 1A and 1B simultaneously then fix formatting
df = tabula.read_pdf('sapp.pdf', 
                     pages='10-13', 
                     multiple_tables=False, 
                     lattice=True, 
                     pandas_options={"header": 0})[0]

df.columns = ['Ordinal#', 'Name', 'CID', 'CAS', '% v/v or g/ml', 'Phase', 'Solvent']

# Fix data from Table 1A
tmp1 = df.iloc[:84].copy().set_index('Ordinal#')

# Add two molecules that got cut-off in PDF read
tmp1.loc['59'] = ['ethyl butyrate', 7762, '105-54-4', 0.01, 'L', 'mineral oil']
tmp1.loc['132'] = ['thiophene', 8030, '110-02-1', 2.5, 'L', 'mineral oil']
tmp1.reset_index(inplace=True)

# Fix data from Table 1B
tmp2 = df.iloc[85:].copy().shift(periods=1, axis='columns')

# Add 1 molecules that got cut-off in PDF read
tmp2.loc[len(tmp2)] = [np.nan, 'ethyl decanoate', 8048, '110-38-3', 100, 'L', np.nan]

odorants = pd.concat([tmp1, tmp2], axis=0)
odorants.CID = odorants.CID.astype(int)
odorants = odorants.set_index('CID').sort_index()
odorants.Phase = odorants.Phase.map({'L': 'liquid', 'S': 'solid'})

odorants.head()

# Get molecule info using CIDs
molecules = pd.DataFrame(pyrfume.from_cids(odorants.index.to_list())).set_index('CID').sort_index()
molecules.head()

# Load perfumer-provided names for each mixture (Supp. Table 2)
perfumer_names = tabula.read_pdf('sapp.pdf', 
                     pages='14', 
                     multiple_tables=False, 
                     lattice=True, 
                     pandas_options={"header": 0})[0]

perfumer_names.iloc[-1:] = perfumer_names.iloc[-1:].shift(periods=1, axis=1)
perfumer_names.ffill(inplace=True)
perfumer_names.columns = ['No. Comp. in Mixture', 'ID task', 'Version 1', 'Version 2', 'Version 3', 'Version 4']
perfumer_names['ID task'] = perfumer_names['ID task'].str.replace('- ', '-')
perfumer_names['Stimulus'] = 'white'

perfumer_names = pd.melt(perfumer_names, id_vars=['Stimulus', 'No. Comp. in Mixture', 'ID task'],
                         var_name='Version', value_name='Descriptors' )

perfumer_names.Version = perfumer_names.Version.str[-1]
perfumer_names.Descriptors = perfumer_names.Descriptors.str.replace(' -\r', '-').str.split('\r').apply(';'.join)         .str.replace(' -', '-').str.replace('/ ', '/')
    
perfumer_names = perfumer_names.set_index(['Stimulus', 'No. Comp. in Mixture', 'Version']).sort_index()

perfumer_names.head()

# Load mixture compositions (Supp. Table 7); the mulipage (p 19-20), multi-column table was unreadable by tabula; 
# had to manually parse as SuppTable.csv to get started
mixtures = pd.read_csv('SuppTable7.csv', header=0)
mixtures = pd.melt(mixtures, var_name='Set', value_name='tmp').dropna()

mixtures.Set = mixtures.Set.str[-1]
mixtures.Set = mixtures.Set.astype(int)
mixtures['Ordinal#'] = mixtures.tmp.str.replace('^', '', regex=False).str.replace('#', '', regex=False)     .str.replace('*', '', regex=False).str.replace('$', '', regex=False).str.replace('&', '', regex=False).str.strip()

# Add CID's by matching to Ordinal#
ord_to_cid = dict(zip(odorants.loc[~odorants['Ordinal#'].isna(), 'Ordinal#'], odorants[~odorants['Ordinal#'].isna()].index))
mixtures['CID'] = mixtures['Ordinal#'].str.strip().map(ord_to_cid)

# Map superscripts to # of components in mixtures
# From Table 7 caption: #: 30, ^: 20, $: 15; *: 10, &: 4, Bold font: 1
superscript_dict = {'#': 30, '^': 20, '$': 15, '*': 10, '&': 4}

for k, v in superscript_dict.items():
    mixtures[str(v) + '-comp'] = mixtures.tmp.apply(lambda x: 1 if k in x else 0)

# Indication for mono-molecule (1-component) if no superscripts or (add manually) bold in the pdf table
mixtures['1-comp'] = mixtures.tmp .apply(lambda x: 0 if ('^' in x) | ('#' in x) | ('$' in x) | ('*' in x) | ('&' in x) else 1)

mixtures.loc[(mixtures.Set == 1) & (mixtures['Ordinal#'] == '47'), '1-comp'] = 1
mixtures.loc[(mixtures.Set == 2) & (mixtures['Ordinal#'] == '81'), '1-comp'] = 1
mixtures.loc[(mixtures.Set == 3) & (mixtures['Ordinal#'] == '138'), '1-comp'] = 1
mixtures.loc[(mixtures.Set == 4) & (mixtures['Ordinal#'] == '38'), '1-comp'] = 1
        
mixtures = pd.melt(mixtures, id_vars=['Set', 'tmp', 'Ordinal#', 'CID'], var_name='No. Comp. in Mixture', value_name='isMember')
mixtures['No. Comp. in Mixture'] = mixtures['No. Comp. in Mixture'].str.split('-').str[0].astype(int)

mixtures = mixtures[mixtures.isMember == 1].groupby(by=['Set', 'No. Comp. in Mixture']).agg({'isMember': 'sum', 'CID': list})
mixtures['Stimulus'] = 'white'
mixtures = mixtures.drop('isMember', axis=1).reset_index().set_index(['Stimulus', 'No. Comp. in Mixture', 'Set'])     .sort_index()

mixtures.head()

# Load 'white' descriptors percentage of applicability data (Supp. Table 3)
descriptors = tabula.read_pdf('sapp.pdf', 
                     pages='15-16', 
                     multiple_tables=False, 
                     lattice=True, 
                     pandas_options={"header": None})[0]

descriptors.columns = ['Idx', 'Descriptors', 'P.A.']
descriptors.dropna(subset=['Idx'], inplace=True)
descriptors['Idx'] = descriptors['Idx'].astype(int)
descriptors.set_index('Idx', inplace=True)

# Add 2 rows that got cut-off in PDF read
descriptors.loc[34] = ['Black pepper', 7.1]
descriptors.loc[107] = ['Sauerkraut', 14.9]

descriptors.Descriptors = descriptors.Descriptors.str.split(',').apply(';'.join)
descriptors['Stimulus'] = 'white'
descriptors = descriptors.set_index('Stimulus').sort_index()

descriptors.head()

# Write to disk
# 'parsed' files
odorants.to_csv('Odorants.csv')

# 'processed' files
molecules.to_csv('molecules.csv')
mixtures.to_csv('stimuli.csv')
perfumer_names.to_csv('behavior_1.csv')
descriptors.to_csv('behavior_2.csv')

