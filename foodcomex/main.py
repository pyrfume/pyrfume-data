#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import requests
import os
import bs4
from tqdm import tqdm
import pyrfume

# Scraping functions
BASE_URL = 'https://foodcomex.org/foodcomex_compounds?c=public_id&d=down&page='

def get_soup(num):
    page_id = str(num)
    try:
        resp = requests.get(BASE_URL + page_id) 
        return bs4.BeautifulSoup(resp.content, 'html.parser')
    except:
        print('No soup for: ', page_id)
        return None
    
def get_compounds(soup):
    compounds = []
    for tr in soup.find('tbody').find_all('tr'):
        tmp = [td for td in tr.find_all('td')]
        index = tmp[0].text
        name = tmp[1].find('a').text
        cas = str(tmp[1]).split('<hr/>')[-1].split('<')[0]
        compounds.append([index, name, cas])
    return compounds

# Scrape to get compounds (Name and CAS) if compounds.csv doesn't exist
if not os.path.isfile('compounds.csv'):
    compounds = []
    for i in tqdm(range(1,50)):
        soup = get_soup(i)
        compounds += get_compounds(soup)

    df = pd.DataFrame(compounds, columns=['Index', 'Name', 'CAS']).set_index('Index')
    df.to_csv('compounds.csv')
else:
    df = pd.read_csv('compounds.csv')

df = df[~df.duplicated()] # Remove duplicates if any
df.CAS = df.CAS.replace('', np.nan)
print(df.shape)
df.head()

# Fetch CIDs using CAS
cas = df.CAS.dropna().to_list()
cids = pyrfume.get_cids(cas)

# Look for missing CIDs by name
df['CID'] = df.CAS.map(cids)
names_missing_cids = df[(df.CID == 0) | (df.CID.isna())]
cids2 = pyrfume.get_cids(names_missing_cids.Name.to_list(), kind='name')

df.loc[names_missing_cids.index, 'CID']  = df.loc[names_missing_cids.index, 'Name'].map(cids2)
print(df[df.CID == 0].shape[0], 'CIDs still missing')

# Get molecule info from CIDs
all_cids = df[df.CID > 0].CID.astype(int).to_list()
molecules = pd.DataFrame(pyrfume.from_cids(all_cids)).set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()] # Remove duplicates
molecules.head()

# Write to disk
molecules.to_csv('molecules.csv')

