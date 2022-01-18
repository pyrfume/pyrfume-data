# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import pandas as pd
import pyrfume
from pyrfume.odorants import from_cids, get_cids

file_name = '12868_2012_2787_MOESM1_ESM.xlsx'
# -

subjects = pd.read_excel(file_name, sheet_name='demographics').set_index('UID')
subjects.index.name = 'subject'
subjects.to_csv('subjects.csv')

odorants = pd.read_excel(file_name, sheet_name='odours and sequence of stimuli')

replacements = {'propylen glycol': 'propylene glycol',
                '2-decenal': 'CC(=CCCC(=CC=O)C)C', # racemic 19801
                'citral': 'CC(=CCCC(=CC=O)C)C', # racemic 8843
                'isobornyl acetate': 'CC(=O)O[C@@H]1C[C@@H]2CC[C@@]1(C2(C)C)C', # 61061 matches isomeric SMILES on GoodScents
               }
odorants['name'] = odorants['odour'].replace(replacements)

# + tags=[]
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')   
names_to_cids = get_cids(odorants['name'].values)
# -

substances = [key for key, value in names_to_cids.items() if value<=0]
for i, substance in enumerate(substances):
    names_to_cids[substance] = -(i+1)

substances = pd.Series(substances, index=range(-1, -len(substances)-1, -1), name='name').to_frame()
substances.to_csv('substances.csv')
substances.head()

molecules = pd.DataFrame(from_cids(list(names_to_cids.values()))).set_index('CID').sort_index()
molecules.to_csv('molecules.csv')
molecules.head()

mixtures = odorants[['odour', 'concentration', 'dilution', 'solvent']].copy()
mixtures = mixtures[odorants['vial #'].notnull()]
mixtures['CID'] = odorants['name'].apply(names_to_cids.get)
mixtures = mixtures[~mixtures.duplicated()].reset_index(drop=True)
mixtures.to_csv('mixtures.csv')
mixtures.head()

thresholds = pd.read_excel(file_name, sheet_name='thresholds').set_index('UID')
thresholds.columns = thresholds.columns.map(names_to_cids.get)
thresholds.columns.name = 'CID'
thresholds.index.name = 'subject'
thresholds.to_csv('behavior_1.csv')
thresholds.head()

descriptors = pd.read_excel(file_name, sheet_name='descriptors', header=[0, 1], index_col=[0, 1])
descriptors.index.names = ['subject', 'rep']
descriptors.columns.names = ['mixture', 'descriptor']
mixture_names = descriptors.columns.levels[0]
mixture_ids = mixture_names.map(names_to_cids.get)
descriptors.columns = descriptors.columns.set_levels(mixture_ids, level=0)
descriptors.to_csv('behavior_2.csv')
descriptors.head()

# +
intensity = pd.read_excel(file_name, sheet_name='intensity and pleasantness',
                          header=[1, 2], index_col=[0, 1]).iloc[:, :138]

def fix_indices(df):
    df.index.names = ['subject', 'rep']
    df.columns.names = ['conc', 'mixture']
    
    concs = df.columns.get_level_values(0)
    odours = df.columns.get_level_values(1)

    # Join duplicate mixtures
    odours = odours.map(lambda x: x.split('.')[0])
    
    # Fix mixture names
    replacements = {'propylen glycol': 'propylene glycol'}  
    odours = odours.map(lambda x: replacements.get(x, x))

    mix_ids = []
    for conc, odour in zip(concs, odours):
        pair = mixtures[mixtures['odour']==odour]
        # Two concentrations of everything except solvent
        assert (pair.shape[0] == 2) or (odour in ['paraffin oil', 'propylene glycol']), \
            print(odour, pair.shape[0])
        conc = conc.replace('medium', 'high')  # Naming covention
        if pair.shape[0]==2:
            mix_id = pair[pair['concentration']==conc].index[0]
        else:
            mix_id = pair.index[0]
        mix_ids.append(mix_id)

    df.columns = pd.Index(mix_ids, name='mixture')
    return df

intensity = fix_indices(intensity)
intensity.to_csv('behavior_3.csv')
intensity.head()
# -

pleasantness = pd.read_excel(file_name, sheet_name='intensity and pleasantness',
                             header=[1, 2], index_col=[0, 1]).iloc[:, 143:]
pleasantness = fix_indices(pleasantness)
pleasantness.to_csv('behavior_4.csv')
pleasantness.head()
