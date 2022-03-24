# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Processing pipeline for Mainland et al, 2015

import pandas as pd
import pyrfume
from pyrfume import get_cids
from pyrfume import from_cids

odorants = pd.read_csv('Odors.tsv', sep='\t')
molecules = pd.DataFrame(from_cids(odorants[~odorants.CID.isna()]['CID'].astype(int))) #Retrieve Pubchem informations from CID

molecules = molecules.merge(odorants, on='CID', how="outer") #Merge on CID, outer merge if mixtures are wanted

molecules = molecules.drop(['SMILES'], axis=1).drop_duplicates(subset=['Odor'])
molecules = molecules.rename(columns={'Odor': 'Odorant'}).set_index("CID") #Change column name for consistency
molecules.to_csv('molecules.csv')

receptors = pd.read_csv('Receptors.tsv', sep='\t')
receptors.set_index('Gene')
receptors.to_csv('receptors.csv')

screening = pd.read_csv('PrimaryScreen.tsv', sep='\t')
screening['OR'] = screening['OR'].fillna(0).astype(int)
screening = screening.rename(columns={'Odor': 'Odorant'})
screening.to_csv('behavior_1.csv')

data = pd.read_excel('sdata20152-s2.xls')

data = data.rename(columns={'Odor': 'Odorant'})
data.to_csv('behavior_2.csv')
