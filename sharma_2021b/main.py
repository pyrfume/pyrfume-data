#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import requests
import bs4
import os
import pyrfume

# Only need to scrape for raw data and concatentate to single .csv file once
# If raw/ doesn't exist; download .xlsx file for each of the 106 primary odors
if not os.path.isdir('raw'):
    os.mkdir('raw')
    
    # scrape mapping for primary odor # to name
    resp = requests.get('https://olfab.iiita.ac.in/olfactionbase/odor') 
    soup = bs4.BeautifulSoup(resp.text, 'html')
    options = soup.find_all('option')
    num_to_primary_odor = {}
    for n in range(106): # 106 primary odors
        num_to_primary_odor[n+1] = options[n].text.replace(' ', '_').replace('/', '')
        
    # Download each .xlsx file to raw/
    for n, name in num_to_primary_odor.items():
        url = "https://olfab.iiita.ac.in/olfactionbase/odors/export?odor=" + str(n) + "&amp;subodor="
        resp = requests.get(url)
        output = open('raw\\' + str(n) + '-' + name + '.xlsx', 'wb')
        output.write(resp.content)
        output.close()

# If odors.csv doesn't exist; create from raw/*.xlsx files
if not os.path.isfile('odors.csv'):
    files = os.listdir('raw')
    frames = []
    for filename in files:
        frames.append(pd.read_excel('raw/' + filename, engine='openpyxl'))

    df2 = pd.concat(frames, axis=0, ignore_index=True)
    df2.index.names = ['Ser. No.']
    df2.rename(columns={'CAS-ID’s': 'CAS Number'}, inplace=True)
    df2['CAS Number'] = df2['CAS Number'].astype(str).str.strip()
    df2.drop('Sr.No.', axis=1).to_csv('odors.csv')
    
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

# Load odorants
odorants = pd.read_csv('odorants.csv', index_col=0)
odorants.rename(columns={'PubChem':'CID'}, inplace=True)
odorants.head()

# Load odorless
odorless = pd.read_csv('odorless.csv', index_col=0)
odorless.rename(columns={'PubChem':'CID'}, inplace=True)
odorless.head()

# Merge CIDs and get molecule info
cids = list(set(odorants['PubChem'].tolist())) + odorless['PubChem'].tolist()
info_dict = pyrfume.from_cids(cids)

# create dataframe for molecules.csv
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.dropna(thresh=2, inplace=True) # One CID returned no info...
molecules.head()

# Load full set of odors and add CIDs
odors = pd.read_csv('odors.csv', index_col=['CAS Number', 'Chemical Name'])
odors.index.names = ['CAS No', 'Name']
odors = odors.join(odorants.set_index(['CAS No', 'Name'])['CID'])
odors.head()

# Create behavior.csv from odor descriptors and odorless chemicals
no_odor = odorless[['Name', 'CAS No', 'CID']].set_index(['CAS No','Name'])
no_odor['Primary Odor'] = 'Odorless'

behavior = pd.concat([odors.drop('Ser. No.', axis=1), no_odor])
behavior['Primary Odor'] = behavior['Primary Odor'].str.replace(' ','')
behavior = behavior.set_index(['CID', 'Primary Odor']).sort_index()
behavior.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')