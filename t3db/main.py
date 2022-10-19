#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import requests
import os
import bs4
import pickle
import pyrfume

# Scraping functions
def get_soup(url):
    try:
        resp = requests.get(url) 
        return bs4.BeautifulSoup(resp.content, 'html.parser')
    except:
        print('Could not get: ', url)
        return None

def parse_soup(soup):
    data = []
    for row in soup.find('tbody').find_all('tr'):
        td = row.find_all('td', limit=2)
        t3d_id = td[0].text
        name = td[1].find('a').text
        cas = str(td[1]).split('<hr/>')[-1].split('<')[0].strip()
        data.append([t3d_id, name, cas])
    return data

base_url = 'http://www.t3db.ca/toxins?c=title&d=up&page='

# Scrape to get list of toxins if scraped_data.csv does not exits
if not os.path.isfile('scraped_data.csv'):
    compounds = []
    for page in range(1,149):
        soup = get_soup(base_url + str(page)) 
        if soup is not None:
            compounds += parse_soup(soup)
        print('Page', str(page), 'done')
        
    df = pd.DataFrame(compounds, columns=['Accession Number', 'Name', 'CAS'])
    df.to_csv('scraped_data.csv')
else:
    df = pd.read_csv('scraped_data.csv')

# In cases where multiple CAS were provided, keep first listed
df.CAS = df.CAS.str.split(' ').str[0].replace('', np.nan)
df.head()

# Fetch CIDs from CAS
cas = df.CAS.dropna().to_list()
cids = pyrfume.get_cids(cas)

# Try to find missing CIDs using names
df['CID'] = df.CAS.map(cids)
names_missing_cids = df[(df.CID == 0) | (df.CID.isna())]
cids2 = pyrfume.get_cids(names_missing_cids.Name.to_list())

print(len([k for k, v in cids2.items() if not v]), 'CIDs still missing')

df.loc[names_missing_cids.index, 'CID'] = df.loc[names_missing_cids.index, 'Name'].map(cids2)

# Get molecule info
all_cids = df.CID.dropna().astype(int).to_list()
molecules = pd.DataFrame(pyrfume.from_cids(all_cids)).set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()] # Remove duplicates
molecules.head()

# Write to disk
molecules.to_csv('molecules.csv')

