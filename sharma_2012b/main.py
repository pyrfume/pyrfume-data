#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import requests
import bs4
import os
import shutil
from time import sleep
import pickle
import sys
import pyrfume

# Only need to scrape for raw data once
# If odors.csv doesn't exist, download .xlsx file for each of the 106 primary odors to create
if not os.path.isfile('odors.csv'):
    # If temp/ doesn't exist, create and download .xlsx files for each primary odor
    if not os.path.isdir('temp'):
        os.mkdir('temp')

        # scrape mapping for primary odor # to name
        resp = requests.get('https://olfab.iiita.ac.in/olfactionbase/odor') 
        soup = bs4.BeautifulSoup(resp.text, 'html')
        options = soup.find_all('option')
        num_to_primary_odor = {}
        for n in range(106): # 106 primary odors
            num_to_primary_odor[n+1] = options[n].text.replace(' ', '_').replace('/', '')

        # Download each .xlsx file to temp/
        for n, name in num_to_primary_odor.items():
            url = 'https://olfab.iiita.ac.in/olfactionbase/odors/export?odor=' + str(n) + '&amp;subodor='
            resp = requests.get(url)
            output = open('temp\\' + str(n) + '-' + name + '.xlsx', 'wb')
            output.write(resp.content)
            output.close()

    # Read in data for each primary odor, concatentate into singel dataframe, and save as odors.csv
    files = os.listdir('temp')
    frames = []
    for filename in files:
        frames.append(pd.read_excel('temp\\' + filename, engine='openpyxl'))

    df = pd.concat(frames, axis=0, ignore_index=True)
    df.index.names = ['Ser. No.']
    df.rename(columns={'CAS-ID’s': 'CAS Number'}, inplace=True)
    df['CAS Number'] = df['CAS Number'].astype(str).str.strip()
    df.drop('Sr.No.', axis=1).to_csv('odors.csv')
    
    # Delete temporary files
    shutil.rmtree('temp')
    
# If odorants.csv doesn't exist, scrape to create
if not os.path.isfile('odorants.csv'):
    odorants = {}
    for page in range(1,161):
        url = 'https://olfab.iiita.ac.in/olfactionbase/chemicals?page=' + str(page)
        resp = requests.get(url) 
        soup = bs4.BeautifulSoup(resp.content, 'html.parser')

        if page == 1:
            col_names = [x.text.strip() for x in soup.find_all('th', scope='col')]

        table = soup.find('tbody')
        rows = table.find_all('tr')
        for row in rows:
            cols = [x.text.strip() for x in row.find_all('td')]
            odorants[cols[0]] = cols[1:]
            
    odorants_df = pd.DataFrame.from_dict(odorants, orient='index', columns=col_names[1:])
    odorants_df.index.names = [col_names[0]]
    odorants_df.drop('Details', axis=1).to_csv('odorants.csv')
    
# If odorless.csv doesn't exist, scrape to create
if not os.path.isfile('odorless.csv'):
    odorless = {}
    for page in range(1,46):
        url = 'https://olfab.iiita.ac.in/olfactionbase/chemicals/odorless?page=' + str(page)
        resp = requests.get(url) 
        soup = bs4.BeautifulSoup(resp.content, 'html.parser')

        # Get column names
        if page == 1:
            col_names = [x.text.strip() for x in soup.find_all('th', scope='col')]

        table = soup.find('tbody')
        rows = table.find_all('tr')
        for row in rows:
            cols = [x.text.strip() for x in row.find_all('td')]
            odorless[cols[0]] = cols[1:]
            
    odorless_df = pd.DataFrame.from_dict(odorless, orient='index', columns=col_names[1:])
    odorless_df.index.names = [col_names[0]]
    odorless_df.drop('Details', axis=1).to_csv('odorless.csv')
    
# If OR-odorant_pairs.csv doesn't exist, scrape to create
if not os.path.isfile('OR-odorant_pairs.csv'):
    # Download .xlsx file of OR-odorant pairs
    url = 'https://olfab.iiita.ac.in/olfactionbase/or-odorant/export?page=1'
    resp = requests.get(url)  
    output = open('temp.xlsx', 'wb')
    output.write(resp.content)
    output.close()
    
    # Read in data from temporary .xlsx file and save as OR-odorant_pairs.csv
    df = pd.read_excel('temp.xlsx', engine='openpyxl')
    df.set_index('Sr. No').to_csv('OR-odorant_pairs.csv')
    
    # Delete temporary files
    os.remove('temp.xlsx')

# Functions to scrape for odor strength, physicochemical properties, and pharmacokinetic profile
def get_soup(sn):
    url = 'https://olfab.iiita.ac.in/olfactionbase/chemical/' + str(sn)
    resp = requests.get(url) 
    return bs4.BeautifulSoup(resp.content, 'lxml')

def get_odor_strength(soup):
    strength = None
    h4 = soup.select_one('h4:-soup-contains("Odor Profile")')
    if h4:
        table = h4.find_next('table')
        td = table.select_one('td:-soup-contains("Strength")')
        if td:
            strength = td.find_next('td').text
    return strength

def get_physicochemical(soup):
    h4 = soup.select_one('h4:-soup-contains("Physicochemical properties")')
    table = h4.find_next('tbody')
    return {tr.th.text: tr.td.text for tr in table.find_all('tr')}

def get_pharmacokinetic(soup):
    h4 = soup.select_one('h4:-soup-contains("Pharmacokinetic profile")')
    table = h4.find_next('tbody')
    return {tr.th.text: tr.td.text for tr in table.find_all('tr')}

# Load odorants
odorants = pd.read_csv('odorants.csv', index_col=0)
odorants.rename(columns={'PubChem':'CID'}, inplace=True)
odorants.head()

# Load odorless
odorless = pd.read_csv('odorless.csv', index_col=0)
odorless.rename(columns={'PubChem':'CID'}, inplace=True)

# Isopulegol (CID = 24585) is listed as both odorant and odorless, drop from odorless
odorless.drop(odorless.loc[odorless['Name']=='Isopulegol'].index, inplace=True)
odorless.head()

# Merge CIDs and get molecule info
cids = list(set(odorants['CID'].tolist())) + odorless['CID'].tolist()
info_dict = pyrfume.from_cids(cids)

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.dropna(thresh=2, inplace=True) # One CID returned no info...
molecules.head()

# If necessary, scrape for odor strength, physicochemical properties, and pharmacokinetic profile
if not os.path.isfile('soups.pkl'):
    # Serial # dictionary for odorants and odorless
    sn_offset = len(odorants) # chem info pages odorless use ser. no. sequentially after last odorant
    sn_dict = dict(zip(odorants.index, odorants['CID']))
    sn_dict.update(dict(zip(odorless.index + sn_offset, odorless['CID'])))

    soups = {}
    # This loop to scrape takes ~5 hrs to complete; pickle soups after to avoid running more than once
    for sn, cid in sn_dict.items():
        # Chlorimide (CID = 76939) has broken link to Details page, have to skip
        if cid not in soups and cid != 76939:
            soups[cid] = get_soup(sn)
            sleep(2) # add delay to avoid server rejecing requests
    
    sys.setrecursionlimit(10000)
    with open('soups.pkl', 'wb') as file:
        pickle.dump(soups, file)

else:    
    with open('soups.pkl', "rb") as file:
        soups = pickle.load(file)

# Load full set of odors and add CIDs
odors = pd.read_csv('odors.csv', index_col=['CAS Number', 'Chemical Name'])
odors.index.names = ['CAS No', 'Name']
odors = odors.join(odorants.set_index(['CAS No', 'Name'])['CID'])
odors.head()

# Create behavior_1.csv from odor descriptors and odorless chemicals
no_odor = odorless[['Name', 'CAS No', 'CID']].set_index(['CAS No','Name'])
no_odor['Primary Odor'] = 'Odorless'
no_odor.head()

# Get odor strength from soups
odor_strength = {}
for cid, soup in soups.items():
    odor_strength[cid] = get_odor_strength(soup)

behav1 = pd.concat([odors.drop('Ser. No.', axis=1), no_odor])
behav1['Primary Odor'] = behav1['Primary Odor'].str.replace(' ','')

# Add odor strength where available
behav1['Strength'] = behav1['CID'].map(odor_strength)
behav1 = behav1.set_index(['CID', 'Primary Odor']).sort_index()
behav1.head()

# Load OR-odorant pairs -> behavior_2.csv
pairs = pd.read_csv('OR-odorant_pairs.csv', index_col='Sr. No')
pairs.rename(columns={'PubChem':'CID', 'Receptor':'Olfactory Receptor'}, inplace=True)
pairs.head()

behav2 = pairs[['CID','Olfactory Receptor']]
behav2 = behav2.set_index('CID').sort_index()
behav2.head()

# Create physics.csv from physicochemical and pharmacokinetic properties
physics_dict = {}
for cid, soup in soups.items():
    physics_dict[cid] = get_physicochemical(soup)
    physics_dict[cid].update(get_pharmacokinetic(soup))

physics = pd.DataFrame.from_dict(physics_dict, orient='index')
physics.index.names = ['CID']
physics.sort_index(inplace=True)

# Drop unneeded columns
drop_list = ['Molecular Weight (g/mol)', 'Mass (g/mol)']
physics.drop(axis=1, columns=drop_list, inplace=True)
physics.head()

# Write to disk
molecules.to_csv('molecules.csv')
behav1.to_csv('behavior_1.csv')
behav2.to_csv('behavior_2.csv')
physics.to_csv('physics.csv')