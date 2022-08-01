#!/usr/bin/env python
# coding: utf-8

import requests
import numpy as np
import urllib.request
import pandas as pd
import bs4
import os
import pyrfume
import pickle

# Scraping functions
def get_soup(url):
    response = requests.get(url)
    return bs4.BeautifulSoup(response.text, 'html')

# Gets list of compounds on each listing page
def get_compounds(soup):
    compounds = []
    rows = soup.find('tbody').find_all('tr')
    for row in rows:
        compounds.append([td.text for td in row.select('td', limit=3)])
    return compounds

# Gets properties for each compound from its page specefied by the FooDB ID
def get_compound_properties(foodb_id):
    try: # Will fail if bad url or no 'Predicted Properties' for this compound
        fp = urllib.request.urlopen(base_url + '/' + foodb_id)
        response = fp.read().decode('utf8')
        fp.close()

        soup = bs4.BeautifulSoup(response, 'html')

        table = soup.find('th',text='Predicted Properties').findNext('tbody')
        data = [[td.findChildren(text=True)[0] for td in tr.find_all('td', limit=2)] for tr in table.find_all('tr')]
        for prop in ['IUPAC name', 'InChI Identifier', 'InChI Key']:
            data.append([prop, soup.find('th',text=prop).findNext('td').text]) 
        message = 'ok'
    except:
        data = []
        message = 'fail'
    
    data.append(['Scrape', message])
    print(foodb_id, message)
    return {row[0]: row[1] for row in data}

# Scrape for list of compounds including FooDB IDs which will index the properties page for each compound; takes ~2.5 hrs
base_url = 'https://foodb.ca/compounds'

if not os.path.isfile('compounds.csv'):
    compounds = []
    for page in range(1,2839): # Pages run from 1 to 2838
        soup = get_soup(base_url + '?page=' + str(page))
        compounds += get_compounds(soup)
    
    df = pd.DataFrame(compounds, columns=['FooDB ID', 'Name', 'CAS'])
    df.to_csv('compounds.csv') # Save after 1st scrape to avoid having to repeat
else:
    df = pd.read_csv('compounds.csv', index_col=0)
    
df.head()

# Scrape the page for properties of each compound; pages are indexed by the FooDB ID's; takes ~48 hrs for full scrape
if not os.path.isfile('properties.csv'):
    prop_dict1 = {}
    for i, foodb_id in enumerate(df['FooDB ID']):
        prop_dict1[foodb_id] = get_compound_properties(foodb_id) 
        if not np.mod(i, 100): # Write to disk every 100 compounds in case of crash
            pd.DataFrame.from_dict(prop_dict1, orient='index').to_csv('properties.csv')
    properties = pd.DataFrame.from_dict(prop_dict1, orient='index')
    properties.to_csv('properties.csv') # Save after 1st scrape to avoid having to repeat
else:
    properties = pd.read_csv('properties.csv', index_col=0)

# Need to continue scraping for any missing IDs; server seems to stop responding after a few ~10k requests...
# Retry any failed attempts and IDs yet to be scraped:

if properties.shape[0] != len(df['FooDB ID']):
    missing_ids = list(set(df['FooDB ID']).difference(properties.index.to_list()))

    if not os.path.isfile('properties2.csv'):
        prop_dict2 = {}
        for i, foodb_id in enumerate(missing_ids):
            prop_dict2[foodb_id] = get_compound_properties(foodb_id) 
            if not np.mod(i, 100): # Write to disk every 100 compounds in case of crash
                pd.DataFrame.from_dict(prop_dict2, orient='index').to_csv('properties2.csv')
        temp = pd.DataFrame.from_dict(prop_dict2, orient='index')
        temp.to_csv('properties2.csv') # Save after 1st scrape to avoid having to repeat
    else:
        temp = pd.read_csv('properties2.csv', index_col=0)

    properties = pd.concat([properties, temp])
    properties.to_csv('properties.csv')
    os.remove('properties2.csv')
else:
    print('Scraping complete.')

# 71 compounds did not have 'Predicted Properites' listed in the database; will drop
properties = properties[properties.Scrape == 'ok']
properties.drop('Scrape', axis=1, inplace=True)
properties.head()

# Get CIDs using name
if not os.path.isfile('cids.pkl'):
    cids = pyrfume.get_cids(df['Name'].to_list())
    with open('cids.pkl', 'wb') as f:
        pickle.dump(cids, f)
else:
    with open('cids.pkl', 'rb') as f:
        cids = pickle.load(f)

# Try to get missing CIDs from CAS
if not os.path.isfile('cids2.pkl'):
    missing_names = [k for k,v in cids.items() if not v]
    missing_cas = df.dropna(subset='CAS', axis=0).loc[df.Name.isin(missing_names)].CAS.to_list()
    
    cids2 = pyrfume.get_cids(missing_cas)
    
    with open('cids2.pkl', 'wb') as f:
        pickle.dump(cids2, f)
else:
    with open('cids2.pkl', 'rb') as f:
        cids2 = pickle.load(f)

# Try to get remainder of missing CIDs using InChI Keys
if not os.path.isfile('cids3.pkl'):
    mask = df[(df['FooDB ID'].isin(properties.index)) & 
              ((df.Name.isin([k for k, v in cids.items() if not v])) | 
              (df.CAS.isin([k for k, v in cids2.items() if not v])))]['FooDB ID'].to_list()

    inChIKeys = properties.loc[mask]['InChI Key'].to_list()

    cids3 = pyrfume.get_cids(inChIKeys)
    
    with open('cids3.pkl', 'wb') as f:
        pickle.dump(cids3, f)
else:
    with open('cids3.pkl', 'rb') as f:
        cids3 = pickle.load(f)

# Get info from cids
all_cids = [v for v in cids.values() if v] + [v for v in cids2.values() if v] + [v for v in cids3.values() if v]

print(df.shape[0] - len(set(all_cids)), 'missing after searching by Name, CAS, and InChI Key')

mol_dict = pyrfume.from_cids(all_cids)

# Create physics.csv
df['CID'] = df.Name.map(cids)
df.loc[df.CID == 0, 'CID'] = df.loc[df.CID == 0, 'CAS'].map(cids2)

id_to_cid = dict(zip(df['FooDB ID'], df.CID))

physics = properties.copy()
physics['CID'] = physics.index.map(id_to_cid)
physics['CID'].fillna(0, inplace=True)
physics.loc[physics.CID == 0, 'CID'] = physics.loc[physics.CID == 0, 'InChI Key'].map(cids3)

# Remove any remaining missing CIDs then drop any duplicates
physics = physics[~(physics.CID == 0)]
physics = physics.set_index('CID').sort_index()
physics.index = physics.index.astype(int)
physics = physics[~physics.index.duplicated()]
physics.head()

# Create dataframe for molecules.csv
molecules = pd.DataFrame(mol_dict).set_index('CID').sort_index()

# Only keep molecules that ended up in physics dataframe and drop duplicates
molecules = molecules.loc[physics.index]
molecules = molecules[~molecules.index.duplicated()]

molecules.head()

print(df.shape[0] - molecules.shape[0], 'orphan molecules after searching for CIDs & dropping duplicates')

# Write to disk
molecules.to_csv('molecules.csv')
physics.to_csv('physics.csv')