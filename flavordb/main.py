#!/usr/bin/env python
# coding: utf-8

import requests
import numpy as np
import pandas as pd
import bs4
import os
import pyrfume
import re
from spellchecker import SpellChecker

# Scraping functions
def get_soup(url):
    try:
        response = requests.get(url)
        return bs4.BeautifulSoup(response.text, 'html')
    except:
        print('GET fail')
        return []

def soup_to_cids(soup):
    return [row.text.strip().split('\n')[1] for row in soup.find('tbody').find_all('tr')]

def soup_to_details(soup):
    details = ['Flavor Profile: ', 'FEMA Flavor Profile: ', 'Taste: ', 'Odor: ']
    if soup:
        return {d: soup.find('th', text=d).find_next('td').text for d in details}
    else:
        return {d: 'no soup' for d in details}

# Scrape first to get list of CIDs for the set of molecules; takes ~15 min
base_url = 'https://cosylab.iiitd.edu.in/flavordb'

# Non-specific search querry to access all molecules
search_query = '/molecules?common_name=&functional_group=&flavor_profile=&fema_flavor=&molecular_weight_from=&h_bond_donors=&h_bond_acceptors=&type=&smile=&page='

if not os.path.isfile('cids.csv'):
    cids = []
    for page in range(1,513): # 512 pages listing molecules in the database
        soup = get_soup(base_url + search_query + str(page))
        cids += soup_to_cids(soup)

    pd.DataFrame(cids, columns=['CID']).to_csv('cids.csv', index=None)
else:
    cids = pd.read_csv('cids.csv', usecols=['CID']).CID.to_list()

# Scrape each molecules details page; pages are indexed by CID; takes ~9 hours
if not os.path.isfile('details.csv'):
    data_dict = {}
    for i, cid in enumerate(cids):
        soup = get_soup('https://cosylab.iiitd.edu.in/flavordb//molecules_details?id=' + cid)
        data_dict[cid] = soup_to_details(soup)
        print(i, 'of', len(cids))
        
    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df.index.name = 'CID'
    df.to_csv('details.csv')
else:
    df = pd.read_csv('details.csv', index_col=0)

# Tidy columns names
df.columns = [col.replace(':', '').strip() for col in df.columns]

def filter_list(my_list):
    char_list = [r'[\[\]"&<>.(:%/)]', '-like']
    tmp = []
    for term in my_list:
        term = re.sub('|'.join(char_list), '', term)
        term = term.strip('-').replace('--', '')
        term = re.sub('\d+', '', term).replace('-', ' ').replace('_', ' ')
        for t in term.split():
            if term:
                tmp.append(t)
    return tmp

# Split odor string into lists of terms
df['Odor'] = df.Odor.fillna('').str.replace('\r\n', ' ').str.replace(',', '').str.split().to_list()
df.Odor = df.Odor.apply(lambda x: list(set(x)))
df.Odor = df.Odor.apply(lambda x: filter_list(x))

# Merge flavor descriptors
df['Flavor'] = df[['Flavor Profile', 'FEMA Flavor Profile']].fillna('').agg(' '.join, axis=1).str.lower()
df.Flavor = df.Flavor.str.replace(',', '').str.split().to_list()
df.Flavor = df.Flavor.apply(lambda x: list(set(x)))
df.Flavor = df.Flavor.apply(lambda x: filter_list(x))

df.drop(['Flavor Profile', 'FEMA Flavor Profile', 'Taste'], axis=1, inplace=True)

# Drop molecules that had no data provided
df = df[~df.isna().all(axis=1)]

print(df.shape)
df.head()

# Check spelling and standardization in descriptor columns

# Build remove list by inspecting output of SpellChecker
remove_words = ['soln', 'iaa', 'mgcu', 'ppb', 'ppm', 'odor', 'p', 'concn', 'm', 'dl', 'ordor', 'peely', 'odored',
                'spectroscopically', 'when', 'the', 'and','andor', 'a', 'yet', 'with', 'without', 'as', 'at', 'any', 'an',
                'almost', 'also', 'are', 'be', 'because', 'becoming', 'but', 'available', 'between', 'acquires', 'aroma',
                'available', 'changing', 'character', 'characterized', 'characteristically', 'commercial', 'concentrated',
                'concentration', 'concentrations', 'contains', 'cut', 'description', 'detection', 'develop', 'develops',
                'especially', 'frequently', 'from', 'good', 'grade', 'grades', 'has', 'have', 'having', 'if', 'in', 'it',
                'less', 'like', 'lies', 'levels', 'may', 'mixture', 'more', 'nearly', 'no', 'non', 'not', 'note', 'of',
                'or', 'other', 'otherwise', 'over', 'peel', 'peels', 'plus', 'pure', 'purified', 'quite', 'remotely',
                'resembles', 'resembling', 'room', 'scent', 'similar', 'smell', 'so', 'some', 'sometimes', 'somewhat',
                'specific', 'suggestive', 'temp', 'that', 'thin', 'threshold', 'times', 'to', 'tone', 'tones', 'traces',
                'type', 'typical', 'upon', 'usually', 'vary', 'very', 'weight', 'aging', 'air', 'approaching',
                'approximately', 'high', 'higher', 'is', 'impart', 'impure', 'impurities', 'impurity', 'finer', 'exhibits',
                'exposure', 'french', 'free', 'liquid', 'midway', 'molecular', 'polish', 'perception', 'plates', 'present',
                'preserve', 'quality', 'rather', 'recalling', 'recognition', 'referred', 'relatively', 'resemble',
                'residual', 'rubbing', 'should', 'stench', 'taste', 'temperature', 'tenacity', 'unique', 'vapors',
                'bruised', 'characteristics', 'differing', 'discernable', 'dilute', 'diluted', 'dilution', 'foods', 'none',
                'occasional', 'on', 'practically', 'reminiscent' 'solution', 'solutions', 'valley', 'diffusive',
                'mouthfeel', 'absolute', 'alter', 'alters', 'arise', 'ball', 'box', 'bios', 'bud', 'certain', 'cloth',
                'cothes', 'cloths', 'coop', 'de', 'does', 'help', 'imparts', 'ingredient', 'notes', 'others', 'several',
                'stable', 'stabilize', 'tastes', 'frutti', 'various', 'within', 'flavor', 'flavoring', 'flavors', 'flavor',
                'forms', 'formulation', 'formulations'
]

# Standardize spelling
spell_repl = {'solventy': 'solvent', 'benzenaldehyde': 'benzaldehyde', 'terebinthinate': 'turpentine',
              'sulfidy': 'sulfide', 'ammoniacal': 'ammonia', 'roselike': 'rose', 'blossomy': 'blossom',
              'methylacetopheneone': 'methylacetophenone', 'trichloromethane': 'chloroform', 'cloves': 'clove',
              'characteritic': 'characteristic', 'peppermintlike': 'peppermint', 'spicey': 'spicy',
              'nuancescitrus': 'citrus', 'diacytl': 'diacetyl', 'rancidity': 'rancid', 'descript': 'nondescript',
              'wood': 'woody', 'alcoholic': 'alcohol', 'aldehydes': 'aldehyde', 'apples': 'apple', 'fat': 'fatty',
              'camphoraceous': 'camphor', 'caramellic': 'caramel', 'butter': 'buttery', 'burning': 'burnt',
              'faintly': 'faint', 'feces': 'fecal', 'fish': 'fishy', 'flower': 'flowery', 'flowers': 'flowery',
              'fruit': 'fruity', 'grass': 'grassy', 'herbaceous': 'herbal', 'hyacinths': 'hyacinth',
              'intensely': 'intense', 'jasmin': 'jasmine', 'mildly': 'mild', 'mint': 'minty', 'moderately': 'moderate',
              'musky': 'musk', 'pears': 'pear', 'phenol': 'phenolic', 'piney': 'pine', 'roses': 'rose', 'rosy': 'rose',
              'sickeningly': 'sickening', 'slightly': 'slight', 'spiced': 'spicy', 'spice': 'spicy', 'sulfurous': 'sulfur',
              'sweetish': 'sweet', 'sweetness': 'sweet', 'tar': 'tarry', 'thymol': 'thyme', 'violets': 'violet',
              'weaker': 'weak', 'winey': 'wine', 'vinous': 'wine', 'drying': 'dry', 'freshly': 'fresh',
              'freshness': 'fresh', 'nail': 'nail polish', 'distinct': 'distinctive', 'oily': 'oil', 'strongly': 'strong',
              'terpentine': 'turpentine', 'logenberry': 'loganberry', 'fleafy': 'leafy', 'turnup': 'turnip',
              'gravey': 'gravy', 'bready': 'bread', 'coumarinic': 'coumarin', 'fruitiness': 'fruity', 'cuminseed': 'cumin',
              'sulfury': 'sulfur', 'aldehyde': 'aldehydic', 'almonds': 'almond', 'beans': 'beany', 'bean': 'beany',
              'beefy': 'beef', 'sweetbitter': 'bittersweet', 'brothy': 'broth', 'camphoreous': 'camphor', 'catty': 'cat',
              'cheesy': 'cheese', 'cherries': 'cherry', 'cinnamic': 'cinnamon', 'citral': 'citrus', 'citric': 'citrus',
              'cream':'creamy', 'decomposing': 'decayed', 'dust': 'dusty', 'eart': 'earthy', 'earth': 'earthy',
              'estery': 'ester', 'fragant': 'fragrant', 'greenery': 'green', 'hawthorne': 'hawthorn', 'herb': 'herbal',
              'hots': 'hot', 'jammy': 'jam', 'juicy': 'juice', 'ketone': 'ketonic', 'laundry': 'laundered',
              'leathery': 'leather', 'malty': 'malt', 'meat': 'meaty', 'medical': 'medicinal', 'medicine': 'medicinal',
              'mentholic': 'menthol', 'metal': 'metallic', 'milk': 'milky', 'moldy': 'mold', 'mossy': 'moss',
              'mothballs': 'mothball', 'mouldy': 'mold', 'must': 'musty', 'naphtha': 'naphthalic', 'nut': 'nutty',
              'nuts': 'nutty', 'oxidation': 'oxide', 'oxidized': 'oxide', 'painty': 'paint', 'peanuts': 'peanut',
              'peppery': 'pepper', 'potatoes': 'potato', 'resinous': 'resin', 'roast': 'roasted', 'root': 'rooty',
              'rotting': 'rotten', 'rubbery': 'rubber', 'rummy': 'rum', 'sassafrass': 'sassafras', 'seed': 'seedy',
              'sick': 'sickening', 'skunky': 'skunk', 'smoke': 'smoky', 'smoked': 'smoky', 'soap': 'soapy',
              'soupy': 'soup', 'sweat': 'sweaty', 'tropica': 'tropical', 'tutti': 'tutti frutti', 'valeric': 'valerian',
              'water': 'watery', 'wax': 'waxy', 'yeasty': 'yeast', 'asprin': 'aspirin'
}

# Valid words not in default dictionary
add_to_dictionary = ['aldehydes', 'terpene', 'quinone', 'blossomy', 'terpenic', 'acetaldehyde', 'camphoraceous',
                     'muguet', 'resinous', 'cananga', 'caramellic', 'ammoniacal', 'acetophenone', 'anethole',
                     'empyreumatic', 'mignonette', 'ketonic', 'thymol', 'cresylic', 'trephy', 'hedonic', 'vanillin',
                     'florentina', 'methylacetophenone', 'nonresidual', 'nerol', 'trimethylamine', 'quinoline', 'isoamyl',
                     'nitrobenzene', 'pyridine', 'mercaptan', 'winey', 'benzaldehyde', 'ylang', 'vinous', 'diacetyl',
                     'nail polish', 'furfural', 'palmarosa', 'thujone', 'castoreum', 'ionone', 'indole', 'opoponax',
                     'neroli', 'valeric', 'tuberose', 'damascone', 'citral', 'ambrette', 'thiamin', 'storax', 'buchu',
                     'moscato', 'thuja', 'diterpene', 'anisic', 'lactone', 'naphthyl', 'ocimene', 'guaiacol', 'linalool',
                     'acrylate', 'alkane', 'naphtha', 'cedarwood', 'cistus', 'propionic', 'oakmoss', 'aldehydic',
                     'labdanum', 'heliotropin', 'naphthelene', 'cinnamate', 'tolu', 'galbanum', 'eucalyptol', 'benzoate',
                     'acetoin', 'alliaceous', 'peptone', 'cedarleaf', 'lovage', 'iactonic', 'borneol', 'vetiver', 'cresol'
]

# modifier list
mods = ['faint', 'heavy', 'mild', 'slight', 'weak', 'penetrating', 'strong', 'intense', 'extremely', 'harsh', 'highly',
        'moderate', 'sharp', 'soft', 'suffocating', 'transient', 'diffuse', 'persistent', 'powerful', 'pronounced',
        'transient'
]

# Run this check for both Odor and Flavor columns
terms = [x for xs in df.Odor for x in xs] + [x for xs in df.Flavor for x in xs]
terms = [t.strip() for t in terms if t not in remove_words]
terms = [spell_repl[t] if t in spell_repl else t for t in terms]

spell = SpellChecker()
spell.word_frequency.load_words(add_to_dictionary)

for w in spell.unknown(terms):
    if w not in spell: print(w, '->', spell.correction(w), '?')

# use .value_counts() to assist in standardizing similar descriptors
pd.DataFrame(terms).value_counts().sort_values().sort_index()

# Final filtering of odor percepts
# functions for filtering of percepts and splitting into decriptors and modifiers 
def to_percepts(odor_list):
    tmp = []
    for term in odor_list:
        # Replace for spelling and standardization
        if term in spell_repl:
            term = spell_repl[term]
        # Only keep if not a modifier or in the to_drop list
        if term not in mods and term not in remove_words:
            tmp.append(term)
    return ';'.join(x for x in list(set(tmp)))

def to_modifiers(odor_list):
    tmp = []
    for term in odor_list:
        if term in mods:
            tmp.append(term)
    return ';'.join(x for x in tmp)    

df['Odor Percepts'] = df.Odor.apply(lambda x: to_percepts(x))
df['Odor Modifiers'] = df.Odor.apply(lambda x: to_modifiers(x))
df['Flavor Percepts'] = df.Flavor.apply(lambda x: to_percepts(x))
df['Flavor Modifiers'] = df.Flavor.apply(lambda x: to_modifiers(x))
df.drop(['Odor', 'Flavor'], axis=1, inplace=True)
df.index.name = 'Stimulus'
df.head()

# Get info from cids
mol_info = pyrfume.from_cids(df.index.to_list())

# Create dataframe for molecules.csv
molecules = pd.DataFrame(mol_info).set_index('CID').sort_index()
molecules.head()

# Check that all molecules are accounted for
assert list(df.index) == list(molecules.index)

# Create dataframe for stimuli.csv; all stimuli are CIDs
stimuli = pd.DataFrame(molecules.index, index=molecules.index, columns=['CID'])
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to disck
molecules.to_csv('molecules.csv')
df.to_csv('behavior.csv')
stimuli.to_csv('stimuli.csv')

