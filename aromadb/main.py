#!/usr/bin/env python
# coding: utf-8

# for full reproducibility, the first step is to run Scrape_AromaDb.ipynb

import pandas as pd
import numpy as np
from spellchecker import SpellChecker
import pyrfume

# after scraping I manually added the "Refined Descriptors" column to handle basic non-uniformity
# of raw descriptors, change English spellings to American, and correct some spelling errors

# read raw data from initial scrape, drop entries with no odor descriptors, drop duplicate rows
dfRaw = pd.read_csv('AromaDb_raw.csv', index_col=0, keep_default_na=False, na_values=[''])     .dropna(subset=['Refined Descriptors']).drop_duplicates(subset='CID')

# Get standard information from PubChem for the CIDs
cids = dfRaw['CID'].tolist()
info_dict = pyrfume.from_cids(cids)

# create dataframe for molecules
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules = molecules.loc[~molecules.index.duplicated()] # Remove duplicate rows
molecules.head()

# spellcheck for descriptors; edit 'Refined Descriptor' column of AromaDb_raw.csv as needed
odors = []
for odor in dfRaw['Refined Descriptors'].tolist():
    temp = odor.split(',')
    for t in temp:
        odors.append(t.strip().lower().replace('-like','').replace(' like',''))

spell = SpellChecker()

for tmp in spell.unknown(odors):
    for w in tmp.split(' '):
        if w not in spell: print(w, '->', spell.correction(w), '?')

# use .value_counts() to assist in standardizing similar descriptors
df_temp = pd.DataFrame(odors)
pd.set_option("display.max_rows", None)
df_temp.value_counts().sort_values().sort_index()

# merge similar descriptors and separate modifies; use results from .value_counts() above

# modifier list
mods = ['faint', 'heavy', 'mild', 'slightly', 'slight', 'weak', 'penetrating']

# replacement dictionary
repl = {'baked potato': 'potato',
        'buttery': 'butter',
        'camphor': 'camphoraceous',
        'caramel': 'caramellic',
        'citrusy': 'citrus',
        'clove': 'cloves',
        'fat': 'fatty',
        'flower': 'floral',
        'flowery': 'floral',
        'fruit': 'fruity',
        'grassy': 'grass',
        'herbal': 'herbaceous',
        'jasmin': 'jasmine',
        'liquorice': 'licorice',
        'minty': 'mint',
        'musk': 'musky',
        'must': 'musty',
        'nut': 'nutty',
        'oil': 'oily',
        'oniony': 'onion',
        'piney': 'pine',
        'rosey': 'rose',
        'rosy': 'rose',
        'smoke': 'smoky',
        'spice': 'spicy',
        'sweetness': 'sweet',
        'vanillic': 'vanilla',
        'vegetative': 'vegatable',
        'wine': 'winey',
        'wood': 'woody'}

# function for final filtering step on odor descriptors
def final_filter(odors, repl_dict, mod_list):
    tmp = []
    for term in odors.split(','):
        term = term.strip().lower().replace('-like','').replace(' like','')
        if term in repl_dict:
            term = repl_dict[term]
        if term not in mod_list:
            tmp.append(term)
    return ';'.join(x for x in tmp)

# function to seprate modifiers into separate column
def separate_modifiers(odors, mod_list):
    tmp = []
    for term in odors.split(','):
        term = term.strip().lower().replace('-like','').replace(' like','')
        if term in mod_list:
            tmp.append(term)
    return ';'.join(x for x in tmp)      

# create dataframe for behavior
behavior = dfRaw.copy().drop(['Molecule Name'], axis=1).set_index('CID').sort_index()
behavior.index.name = 'Stimulus'

# final stage of filtering for descriptors
behavior['Filtered Descriptors'] = behavior.apply(lambda row: final_filter(row['Refined Descriptors'], repl, mods), axis=1)

# separate modifiers into separate column
behavior['Modifiers'] = behavior.apply(lambda row: separate_modifiers(row['Refined Descriptors'], mods), axis=1)

behavior.drop(['Refined Descriptors'], axis=1).head()

behavior[behavior.index.duplicated(keep=False)].shape

# Create dataframe for stimulus.csv; All stimuli are CID
stimuli = pd.DataFrame(molecules.index, index=molecules.index.tolist())
stimuli.index.name = 'Stimulus'
stimuli.head()

# write to disk
molecules.to_csv('molecules.csv')
behavior.drop(['Refined Descriptors'], axis=1).to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')
