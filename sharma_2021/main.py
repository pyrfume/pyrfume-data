#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import tabula
from spellchecker import SpellChecker
import pyrfume

# Read data from .pdf source
df = tabula.read_pdf('ci0c01288_si_001-1.pdf', 
                     pages='2-127', 
                     multiple_tables=False, 
                     lattice=True, 
                     pandas_options={"header": 1})[0].drop('Source Database', axis=1).set_index('Sr. No')

# Replace IUPAC name with common name if IUPAC name is missing or contains molecular weight or other non-name entries
df['IUPAC Name'].fillna(df['Common Name'], inplace=True)
df.loc[df['IUPAC Name'].str.contains('g/mol'), 'IUPAC Name'] = df['Common Name']
df.loc[df['IUPAC Name'].str.contains('IUPAC'), 'IUPAC Name'] = df['Common Name']
df.loc[df['IUPAC Name'] == '0', 'IUPAC Name'] = df['Common Name']
df.loc[df['IUPAC Name'].str.contains('?', regex=False), 'IUPAC Name'] = df['Common Name']

# Remove carriage returns and extraneous spaces
def clean_names(name):
    return name.replace('\r', '').replace(' )', ')').replace(' ,', ',').replace('  ', ' ').replace('- ', '-').replace(' -', '-')

df['SMILES'] = df['SMILES'].str.replace('\r', '').str.replace(' ', '')
df['Common Name'] = df['Common Name'].apply(lambda x: clean_names(x))
df['IUPAC Name'] = df['IUPAC Name'].apply(lambda x: clean_names(x))

# Do similar for percepts while als preserving compound descriptors
def clean_percepts(percepts, compound_terms):
    percepts = percepts.replace('\r', ' ')
    for term in compound_terms:
        if term in percepts:
            percepts = percepts.replace(term, term.replace(' ', '_'))
    return [p.replace('_', ' ') for p in percepts.split()]

compound_terms = ['bois de rose', 'cut grass', 'dark chocolate', 'dead animal', 'deep fried', 'blue cheese', 'rotton eggs',
                  'rotten fish', 'fir cone- like', 'tutti frutti', 'graham cracker', 'new mown hay', 'lily of the valley',
                  'old wood', 'passion fruit', 'green pea', 'peanut butter', 'raw potato', 'root beer', 'rotting fruit', 
                  'sandal wood', 'swimming pool', 'unripe fruit', 'unripe banana', 'violet leaf', 'rotten vegetables']

df['Smell Percepts'] = df['Smell Percepts'].apply(lambda x: clean_percepts(x, compound_terms))

# Merge precepts for duplicate odorants (i.e. same common name, IUPAC name, and SMILES)
df = df.groupby(by=['Common Name', 'IUPAC Name', 'SMILES']).agg(sum)
df.reset_index(inplace=True)

df.head()

# Use IUPAC name to get CIDs
cids1 = pyrfume.get_cids(df['IUPAC Name'].to_list())
df['CID'] = df['IUPAC Name'].map(cids1)

# Use common names to look for missing CIDs
mask = df.CID == 0
cids2 = pyrfume.get_cids(df.loc[mask, 'Common Name'].to_list())
df.loc[mask, 'CID'] = df.loc[mask, 'Common Name'].map(cids2)

# Use SMILES to look for still missing CIDs
mask = df.CID == 0
cids3 = pyrfume.get_cids(df.loc[mask, 'SMILES'].to_list())
df.loc[mask, 'CID'] = df.loc[mask, 'SMILES'].map(cids3)

dups = df[df.duplicated('CID', keep=False)].sort_values('CID').copy()

dups['cid1'] = dups['Common Name'].map(pyrfume.get_cids(dups['Common Name'].to_list()))
dups['cid2'] = dups['IUPAC Name'].map(pyrfume.get_cids(dups['IUPAC Name'].to_list()))
dups['cid3'] = dups['SMILES'].map(pyrfume.get_cids(dups['SMILES'].to_list()))

dups.head()

# Manually looked up CID duplicates on Pubchem to disambiguiate; -1 means that couldn't identifty as unique, will drop
cid_replacements = { # df index: new CID
    1600: 49,
    148: 111044,
    2239: 31278,
    69: 6950274,
    2136: 7650,
    153: 11094875,
    154: 11378500,
    136: 1715134,
    132: 1550260,
    2455: 5370448,
    3361: 356,
    217: 6708682,
    1942: 6999975,
    2168: 12304610,
    594: 5325911,
    3120: 62151,
    2190: 23623651,
    3704: 11746218,
    2169: 1201528,
    3842: 6951711,
    3962: -1,
    3366: 6997335,
    681: 5978897,
    3340: 929336,
    1162: 61696,
    1338: 71587638,
    2247: 1715213,
    2248: 1549660,
    2176: 23623868,
    744: 11196,
    2993: 7137,
    3615: 61129,
    830: -1,
    131: 11789380,
    3386: 5869600,
    2407: 6432008,
    1656: -1,
    2547: -1,
    1060: 6028462,
    1420: 90899,
    2204: 7664,
    1884: -1,
    3885: 61875,
    735: 119847,
    840: 62002,
    676: 118487,
    3017: 15613,
    1874: 576026,
    1466: 6021887,
    1686: 61875,
    134: 6430321,
    393: 205986,
    33: 91354,
    2908: 5352470,
    3977: 6429290,
    3885: -1,
    251: 6436548,
    2947: 62387}

for ind, cid in cid_replacements.items():
    df.loc[ind, 'CID'] = cid
    
# Drop where CID was set to -1 since these could not be identified as unique odorants
mask = df.CID == -1
df.drop(df[mask].index, axis=0, inplace=True)
df.head()

# Get info from CIDs
mol_dict = pyrfume.from_cids(df['CID'].tolist())

# Create dataframe for molecules.csv
molecules = pd.DataFrame(mol_dict).set_index('CID').sort_index()
molecules.head()

# Check percepts for spelling
percepts = df['Smell Percepts'].to_list()
percepts = [x.strip('-') for xs in percepts for x in xs]

spell = SpellChecker()
words_to_add = ['oakmoss', 'cananga', 'indole', 'eugenol', 'linalool', 'elemi', 'pyrazine', 'petitgrain', 'tagette',
               'vetiver', 'chlorinous', 'cinnamyl', 'camphoraceous', 'farnesol', 'estery', 'tiglate', 'caprylic',
               'storax', 'davana', 'labdanum', 'boronia', 'ylang', 'neroli', 'cistus', 'erogenic', 'durain', 'agarwood',
               'starfruit', 'ambrette', 'acetoin', 'quinoline', 'damascone', 'isojasmone', 'pyridine', 'buchu', 'muguet',
               'tuberose', 'galbanum', 'dihydrojasmone', 'tolu', 'filbert', 'pocpcorn', 'linalyl', 'immortelle', 'carvone',
               'thujonic', 'palmarosa', 'borneol', 'lactonic', 'brothy', 'jasmone', 'alliaceous', 'aldehydic', 'ketonic',
               'vinous', 'citral', 'castoreum', 'origanum', 'mentholic', 'coumarinic', 'ammonical', 'costus', 'anisic',
               'mercaptan', 'acrylate', 'guaiacwood', 'guaiacol', 'furfural', 'piperidine', 'benzaldehyde', 'opoponax',
               'naphthyl', 'terpenic', 'resinous', 'lovage', 'peely', 'roasty', 'alkane', 'diffusive', 'bready', 'winey']
spell.word_frequency.load_words(words_to_add)

for tmp in spell.unknown(percepts):
    for w in tmp.split():
        if w not in spell: print(w, '->', spell.correction(w), '?')

spell_replace = {
    'cashewnut': 'cashew',
    'heral': 'herbal',
    'boxtree': 'box tree',
    'juniperberry': 'juniper',
    'longifolone': 'longifolene',
    'tearose': 'tea rose',
    'mangnolia': 'magnolia',
    'coaltar': 'coal tar',
    'beewax': 'beeswax',
    'sassafrass': 'sassafras',
    'brocolli': 'broccoli',
    'nail-polish': 'nail polish',
    'rneedle': 'fir needle',
    'sweetpea': 'sweet pea',
    'lifitng': 'lifting',
    'bellpepper': 'bell pepper',
    'pocpcorn': 'popcorn',
    'ambergis': 'ambergris'}

for w1, w2 in spell_replace.items():
    idx = np.where(np.array(percepts) == w1)[0]
    for i in idx:
        percepts[i] = w2

# use .value_counts() to assist in standardizing similar descriptors
percepts_df = pd.DataFrame(percepts)
percepts_df.value_counts().sort_values().sort_index()

# Merge similar descriptors and separate modifies; use results from .value_counts() above
# Modifier list
modifiers = ['faint', 'intense', 'light', 'mild', 'powerful', 'penetrating',
            'strong', 'weak', 'lifting', 'sharp', 'soft', 'suffocating']

# Drop list for artifacts from pervious manipulations
to_drop = ['fi', 'rindy']

# replacement dictionary
standardize_dict = {
    'bean': 'beany',
    'beef': 'beefy',
    'butter': 'buttery',
    'chlorinous': 'chlorine',
    'cut grass': 'grassy',
    'fat': 'fatty',
    'fish': 'fishy',
    'flower': 'floral',
    'fruit': 'fruity',
    'gree': 'green',
    'herbal': 'herbaceous',
    'leaf': 'leafy',
    'meat': 'meaty',
    'metal': 'metallic',
    'mint': 'minty',
    'nut': 'nutty',
    'oil': 'oily',
    'roast': 'roasted',
    'roasty': 'roasted',
    'rubber': 'rubbery',
    'spice': 'spicy',
    'vinous': 'winey',
    'wood': 'woody',
    'bois de rose': 'rosewood',
    'sandal wood': 'sandalwood',
    'fir cone- like': 'fir cone',
    'peanut buttery': 'peanut butter',
    'raw potato': 'potato'}

# functions for filtering of percepts and splitting into decriptors and modifiers 
def to_descriptors(percepts, spell_replace, standardize_dict, to_drop):
    tmp = []
    for term in list(set(percepts)): # This will also remove duplicate terms
        # Replace spelling
        if term in spell_replace:
            term = spell_replace[term]
        # Replace for other standardizations
        if term in standardize_dict:
            term = standardize_dict[term]
        # Only keep if not a modifier or in the to_drop list
        if term not in modifiers and term not in to_drop:
            tmp.append(term)
    return ';'.join(x for x in tmp)

def to_modifiers(percepts, modifiers):
    tmp = []
    for term in list(set(percepts)):
        if term in modifiers:
            tmp.append(term)
    return ';'.join(x for x in tmp)      

# Convert smell percepts into decriptors + modifiers for behavior.csv
behavior = df[['CID', 'Smell Percepts']].copy().set_index('CID').sort_index()
behavior['Decriptors'] = behavior['Smell Percepts'].apply(lambda x: to_descriptors(x, spell_replace, standardize_dict, to_drop))
behavior['Modifiers'] = behavior['Smell Percepts'].apply(lambda x: to_modifiers(x, modifiers))
behavior.drop('Smell Percepts', axis=1, inplace=True)
behavior.head()

# Write to disk
molecules.to_csv('molecules.csv')
behavior.to_csv('behavior.csv')