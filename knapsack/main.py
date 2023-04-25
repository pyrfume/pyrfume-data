# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
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

import numpy as np
import pandas as pd
import requests
import os
import bs4
from tqdm import tqdm
import pyrfume
from pyrfume.odorants import hash_smiles
from time import sleep


# +
# Functions
def get_soup(num):
    db_id = 'C' + str(num).zfill(8)
    try:
        resp = requests.get(BASE_URL + db_id) 
        return bs4.BeautifulSoup(resp.content, 'html.parser')
    except:
        print('No soup for: ', db_id)
        return None
    
def parse_soup(soup):
    details = ['Name', 'Mw', 'CAS RN', 'InChIKey', 'SMILES']
    values = []
    if 'Input key word error' not in soup.text:
        for d in details:
            if d == 'Name':
                value = str(soup.find('th', text='Name').find_next('td')).split('>')[1].split('<')[0].strip()
            else:
                value = soup.find('th', text=d).find_next('td').text
            values.append(value)
        return dict(zip(details, values))
    else:
        print('Invalid soup')
        return None
    
def chunks(my_list, n):
    for i in range(0, len(my_list), n):
        yield my_list[i:i + n]


# +
# Scrape to get data if scrape_data.csv does not exits
BASE_URL = 'http://www.knapsackfamily.com/knapsack_core/information.php?mode=r&word='

if not os.path.isfile('scrape_data.csv'):
    scrape_dict = {}
    for i in tqdm(range(1,50410)):
        soup = get_soup(i)
        scrape_dict['C' + str(i).zfill(8)] = parse_soup(soup)
        
    df = pd.DataFrame.from_dict(scrape_dict, orient='columns').T
    df.index.name = 'Knapsack ID'
    df.to_csv('scrape_data.csv')
else:
    print('loading from file')
    df = pd.read_csv('scrape_data.csv', index_col=0)

# +
df.replace('', np.nan, inplace=True)
df.rename(columns={'CAS RN': 'CAS'}, inplace=True)
df.index.name = 'Knapsack ID'

# Drop rows which have no name, CAS, or SMILES
df.dropna(subset=['Name', 'CAS', 'SMILES'], how='all', inplace=True)

print(df.shape)
df.head()

# +
# Fetch CIDs
chunk_size = 1000
sleep_time = 30

smi = df[~df.SMILES.isna()].SMILES.to_list()
# pd.DataFrame(smi).to_csv('smiles.csv', index=False, header=False)

# cids1 = {}
# for chunk in tqdm(chunks(smi, chunk_size)):
#     cids1.update(pyrfume.get_cids(chunk, kind='smiles'))
#     sleep(sleep_time)

# +
# The above search by SMILES won't complete w/o getting timeout error from PubChem, so instead save SMILES to .csv file and
# run through the PubChem exchange service (https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi), save the results and
# load here:

cids1 = pd.read_csv('smiles_to_cids.txt', sep='\t', index_col=None, header=None)
cids1 = dict(zip(cids1[0], cids1[1]))

df['CID'] = df.SMILES.map(cids1)
# -

# Now try by CAS
cas_missing_cid = df[(df.CID == 0) | (df.CID.isna())]
cids2 = {}
for chunk in tqdm(chunks(cas_missing_cid.CAS.dropna().to_list(), chunk_size)):
    cids2.update(pyrfume.get_cids(chunk))
    sleep(sleep_time)

# Now try by name
df.loc[cas_missing_cid.index, 'CID'] = df.loc[cas_missing_cid.index, 'CAS'].map(cids2)
names_missing_cid = df[(df.CID == 0) | (df.CID.isna())]
cids3 = {}
for chunk in tqdm(chunks(names_missing_cid.Name.dropna().to_list(), chunk_size)):
    cids3.update(pyrfume.get_cids(chunk), kind='name')
    sleep(sleep_time)

# +
# Check how many are still missing
df.loc[names_missing_cid.index, 'CID'] = df.loc[names_missing_cid.index, 'Name'].map(cids3)
still_missing_cid = df[(df.CID == 0) | (df.CID.isna())]

# Only keep those that have SMILES
still_missing_cid = still_missing_cid[~still_missing_cid.SMILES.isna()]

still_missing_cid.drop(['CAS', 'InChIKey', 'CID'], axis=1, inplace=True)
still_missing_cid.columns = ['name', 'MolecularWeight', 'IsomericSMILES']
still_missing_cid['CID'] = [-i for i in range(1, len(still_missing_cid) + 1)]
still_missing_cid.set_index('CID', inplace=True)

print(still_missing_cid.shape)
still_missing_cid.head()
# -

# Get info for molecules.csv
all_cids = df.loc[(df.CID != 0) & (~df.CID.isna()), 'CID'].astype(int)
all_cids = list(set(all_cids))
molecules_tmp = pd.DataFrame(pyrfume.from_cids(all_cids)).set_index('CID').sort_index()

# +
molecules = pd.concat([molecules_tmp, still_missing_cid], axis=0)

# Replace negative integer CID with hash of SMILES
molecules.index = molecules.apply(lambda row: hash_smiles(row['IsomericSMILES']) if row.name < 0 else row.name, axis=1)

# Remove any duplicates
molecules = molecules[~molecules.index.duplicated()].sort_index()

print(molecules.shape)
molecules.head()
# -

# Write to disk
molecules.to_csv('molecules.csv')
