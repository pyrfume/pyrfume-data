#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import requests
import os
import bs4
import pickle
from spellchecker import SpellChecker
import pyrfume

# Function to scrape each odorant page
def get_odorant(base_url, page):
    resp = requests.get(base_url + page) 
    soup = bs4.BeautifulSoup(resp.content, 'html.parser')
    
    # Odorant names
    names = soup.find('p', {'class': 'title'}).text.replace('\n', ', ').split(', ')
    names = [n for n in names if n != ''] # Drop blank if only one name provided
    
    # Odor descriptors
    odors = soup.select_one('p:-soup-contains("Percepts")').text.replace('Percepts: ', '').split(', ')
    
    # CAS and MW
    tmp = soup.find('p', {'class': 'cas'}).text.split('\xa0')
    cas, mw = tmp[0].split(' ')[1], tmp[-1].split(' ')[1]
    
    return {'CAS': cas, 'MW': mw, 'Odors': odors, 'Names': names}

# Get list of webpages for odorants from the base url
base_url = 'https://www.flavornet.org/info/'
resp = requests.get(base_url) 
soup = bs4.BeautifulSoup(resp.content, 'html.parser')
table = soup.find('table')
pages =  [t.text for t in table.find_all('a') if 'html' in t.text]

# Scrape to get data for all odorants
# Pickle after 1st full scrape to avoid having to repeat
if not os.path.isfile('scraped_data.pkl'):
    data_dict = {}
    for page in pages:
        cas = page.split('.')[0]
        
        try:
            data_dict[cas] = get_odorant(base_url, page)
        except:
            data_dict[cas] = {'CAS': None, 'MW': None, 'Odors': None, 'Names': 'No data returned'}
    
    with open('scraped_data.pkl', 'wb') as file:
        pickle.dump(data_dict, file)

else:    
    with open('scraped_data.pkl', "rb") as file:
        data_dict = pickle.load(file)

# As sanity check, I compared the CAS# from the list of pages to what was actually on the page and all match
df = pd.DataFrame.from_dict(data_dict, orient='index')
df['Primary Name'] = df['Names'].str[0]
# Fix name formatting issue for one odorant
df.loc['27576-03-0', 'Primary Name'] = 'p, a -dimethylstyrol'
df.head()

# Get CIDs from CAS
cas = df['CAS'].tolist()
cids = pyrfume.get_cids(cas)

# Try to find missing CID's by name
missing_cas = [k for k, v in cids.items() if not v]
print(len(missing_cas))

missing_names = df['Primary Name'].loc[df.index.isin(missing_cas)].to_list()
cids2 = pyrfume.get_cids(missing_names)

still_missing = [k for k, v in cids2.items() if not v]
still_missing

# Manually add missing CIDs; search by name on PubChem
man_add = {'p, a -dimethylstyrol': 33936,  
         'dehydrocarveol': 527142,
         'p-mentha-dien-hydroperoxide': 564367,
         'methylethylpyrazine': 26332,
         'isopiperitone': 6987,
         '(R)-linden ether': 21606187,
         'epoxy-2-nonenal': 527305,
         'epoxy-2-undecenal': 6429298}

for k, v in man_add.items():
    cids2[k] = v

# Add CIDs to df to indentify duplicates
df['CID'] = df['CAS'].map(cids)
mask = df.CID == 0
df.loc[mask, 'CID'] = df.loc[mask, 'Primary Name'].map(cids2)

# Find duplicate CIDs
pd.set_option('display.max_colwidth', None)
dup = df[df.duplicated('CID', keep=False)]
dup.sort_values('CID', axis=0).head(dup.shape[0])

# Get info from cids
all_cids = df['CID'].tolist()
mol_dict = pyrfume.from_cids(all_cids)

# Create dataframe for molecules.csv
molecules = pd.DataFrame(mol_dict).set_index('CID').sort_index()
molecules = molecules[~molecules.index.duplicated()] # Remove duplicate rows
molecules.head()

# Check odor descriptors for spelling
odors = df['Odors'].to_list()
odors = [x for xs in odors for x in xs]

spell = SpellChecker()

for tmp in spell.unknown(odors):
    for w in tmp.split(' '):
        if w not in spell: print(w, '->', spell.correction(w), '?')

spell_repl = {'terpentine': 'turpentine'}
for w1, w2 in spell_repl.items():
    idx = np.where(np.array(odors) == w1)[0]
    for i in idx:
        odors[i] = w2

# use .value_counts() to assist in standardizing similar descriptors
df_temp = pd.DataFrame(odors)
df_temp.value_counts().sort_values().sort_index()

# Merge similar descriptors and separate modifies; use results from .value_counts() above
# Replacement dictionary
repl = {'apple. rose': ['apple, rose'],
        'weak spice': 'spice'} # only this one odor had a modifier, so just remove 'weak'

# Function to filter odors
def filter_odors(odors, spell_repl, repl):
    tmp = []
    for term in list(set(odors)):
        # Replace spelling
        if term in spell_repl:
            term = spell_repl[term]
        # Replace for other standardizations
        if term in repl:
            term = repl[term]
        if isinstance(term, list):
            tmp += term
        else:
            tmp.append(term)
    return ';'.join(x for x in tmp)

# Odor descriptors -> behavior.csv
behavior = df[['CID', 'Odors']].copy().set_index('CID').sort_index()

# Merge descriptors for dupliate CIDs
behavior = behavior.groupby(level=0).agg(sum)

# Apply filtering function to odor descriptors
behavior['Descriptors'] = behavior.apply(lambda row: filter_odors(row['Odors'], spell_repl, repl), axis=1)
behavior.drop('Odors', axis=1, inplace=True)
behavior.index.name = 'Stimulus'
behavior.head()

# Create dataframe for stimuli.csv; all simuli are CIDs
stimuli = pd.DataFrame(molecules.index, index=molecules.index)
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')

