# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import pandas as pd
import numpy as np
import pyrfume
import camelot

# Tabular data from Supporting Information
pdf_file = 'pnas.1713026115.sapp52-70.pdf'

# Mapping of in-lab codes to IUPAC names, from Tables S1 and S3
code_to_iupac = {
    '1': '3-Methylcyclopentadecan-1-one',
    '(R)-1': '(R)-3-Methylcyclopentadecan-1-one',
    '2': 'Cyclopentadecanone',
    '3': 'Cyclopentadecanol',
    '4': 'Cyclopentadecane',
    '5': 'Thiacyclopentadecane 1-oxide',
    '6': 'Thiacyclopentadecane',
    '7': 'Cyclopentadecanethiol',
    '8': 'Cyclohexadecanone',
    '9': 'Oxacyclohexadecan-2-one',
    '10': 'Cyclohexadec-5-en-1-one',
    '11': 'Cycloheptadecanone',
    '12': 'Cycloheptadec-9-en-1-one',
    '13': '1-Oxacycloheptadec-7-en-2-one',
    '14': '1,4-Dioxacycloheptadecane-5,17-dione',
    '15': 'Cyclododecanone',
    '16': 'Oxacyclotridecan-2-one',
    '17': 'Thiacyclotridecane-1-oxide',
    '18': 'Azacyclotridecan-2-one',
    '19': '1-(4-tert-Butyl-2,6-dimethyl-3,5-dinitrophenyl)ethanone',
    '20': '1-tert-Butyl-3,5-dimethyl-2,4,6-trinitrobenzene',
    '21': '1-tert-Butyl-3,4,5-trimethyl-2,6-dinitrobenzene',
    '22': '1-tert-Butyl-2-methoxy-4-methyl-3,5-dinitrobenzene',
    '23': '2-Methyl-1,3,5-trinitrobenzene',
    '24': '1-Methyl-2,4-dinitrobenzene',
    '25': '1-(1,1,2,6-Tetramethyl-3-propan-2-yl-2,3-dihydroinden-5-yl)ethanone',
    '26': '4,6,6,7,8,8-Hexamethyl-1,3,4,7-tetrahydrocyclopenta[g]isochromene',
    '27': '1-(3,5,5,6,8,8-Hexamethyl-6,7-dihydronaphthalen-2-yl)ethanone',
    '28': '(R)-6,6-Difluoro-3-methylcyclopentadecan-1-one',
    '29': 'R)-7,7-Difluoro-3-methylcyclopentadecan-1-one',
    '30': '(R)-8,8-Difluoro-3-methylcyclopentadecan-1-one',
    '31': '(R)-10,10-Difluoro-3-methylcyclopentadecan-1-one',
    '32': '(R)-3-Methylcyclopentadec-6-en-1-one E:Z = 10:1',
    '(Z)-32': '(R,Z)-3-Methylcyclopentadec-6-en-1-one',
    '33': '(R)-6,6-Difluoro-3-methylcyclopentadec-8-en-1-one E:Z = 3:2',
    '(E)-34': '(R,E)-7,7-Difluoro-3-methylcyclopentadec-10-en-1-one',
    '(Z)-34': '(R,Z)-7,7-Difluoro-3-methylcyclopentadec-10-en-1-one',
    '35': '(R)-8,8-Difluoro-3-methylcyclopentadec-5-en-1-one E:Z = 10:1',
    '(E)-36': '(R,E)-9,9-Difluoro-3-methylcyclopentadec-6-en-1-one',
    '(E)-37': '(R,E)-10,10-Difluoro-3-methylcyclopentadec-6-en-1-one',
    '(Z)-37': '(R,Z)-10,10-Difluoro-3-methylcyclopentadec-6-en-1-one'
}

# +
# Table S2 EC50 data; NaN's indicate "little or no response where no meaningful EC50 or top values can be obtained"
# OR5AN1 responses normalized to highest value of racemic muscone (code = 1), OR1A1 responses to musk ambrette (code = 22)
s2 = camelot.read_pdf(pdf_file, pages='12', flavor='stream')[0].df
s2 = s2.iloc[5:]
s2.columns = ['Ring', 'Odorant', 'Code', 'OR5AN1', 'OR1A1']
s2 = s2.drop('Ring', axis=1).replace('', np.nan).dropna(subset='Odorant').reset_index(drop=True)

# Split EC50 from "top" values (\n) and comparisons to previously published values {/}
s2['OR5AN1'] = s2['OR5AN1'].dropna().apply(lambda x: str(x).split('\n')[0].split('/')[0])
s2['OR1A1'] = s2['OR1A1'].dropna().apply(lambda x: str(x).split('\n')[0].split('/')[0])

# Add IUPAC names
s2['IUPAC Name'] = s2['Code'].map(code_to_iupac)

s2.head()
# -

names = s2['IUPAC Name'].to_list()
cids = pyrfume.get_cids(names, kind='name')
# Can't find CIDs for Thiacyclopentadecane 1-oxide or Thiacyclotridecane-1-oxide (both synthesized in-house)

# +
# Table S4, EC50 data for OR5AN1 responses to luorinated (R)-muscone–related analogs
# The responses are normalized to the highest value of (R)-muscone
s4 = camelot.read_pdf(pdf_file, pages='15', flavor='stream')[0].df
s4.columns = ['Odorant', 'Code', 'EC50', 'Top']
s4 = s4.iloc[5:].drop('Top', axis=1).reset_index(drop=True)

# Add IUPAC names
s4['IUPAC Name'] = s4['Code'].map(code_to_iupac)

s4.head()

# +
# Merge EC50 data -> behavior.csv
behavior = pd.concat([
    s2.set_index('Code')[['OR5AN1', 'OR1A1']].melt(
        var_name='Subject', value_name='EC50 (uM)', ignore_index=False
    ).reset_index(),
    s4.set_index('Code').rename(columns={'EC50': 'OR5AN1'})['OR5AN1'].to_frame().melt(
        var_name='Subject', value_name='EC50 (uM)', ignore_index=False
    ).reset_index()

])
behavior = behavior.rename(columns={'Code': 'Stimulus'}).set_index(['Stimulus', 'Subject']).sort_index()

behavior.head()
# -

# Get CIDS
names = list(code_to_iupac.values())
cids = pyrfume.get_cids(names, kind='name')

# Manually add missing CIDs
# Can't find CIDs for most compounds synthesized in-house
cids['(R)-3-Methylcyclopentadecan-1-one)'] = 7160826 # (R)-muscone

# Molecule info
molecules = pd.DataFrame(pyrfume.from_cids(list(cids.values()))).set_index('CID').sort_index()
molecules.head()

# +
# Stimuli: use in-lab odorant code as stimulus ID
stimuli = pd.concat([
    s2.set_index('Code')[['Odorant', 'IUPAC Name']].copy(),
    s4.set_index('Code')[['Odorant', 'IUPAC Name']].copy()
]).sort_index()
stimuli['CID'] = stimuli['IUPAC Name'].map(cids)
stimuli.CID.replace(0, None, inplace=True)

stimuli.head()
# -

# Subjects
subjects = pd.DataFrame([['OR5AN1', 'OR5AN1'], ['OR1A1', 'OR1A1']], columns=['Subject', 'OR'])
subjects = subjects.set_index('Subject').sort_index()
subjects.head()

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
subjects.to_csv('subjects.csv')
behavior.to_csv('behavior.csv')
