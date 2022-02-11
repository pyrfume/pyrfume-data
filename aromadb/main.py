#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import pyrfume

# after scraping I manually added the "Refined Descriptors" column to handle non uniformity of raw descriptors
# read raw data from initial scrape, drop entries with no odor descriptors
dfRaw = pd.read_csv('AromaDb_raw.csv', index_col=0).dropna(subset=['Raw Descriptors'])


# Get standard information from PubChem for the CIDs
cids = dfRaw['CID'].tolist()
info_dict = pyrfume.from_cids(cids)

# create dataframe for molecules
molecules = pd.DataFrame(info_dict).set_index('CID').sort_index()
molecules.head()

# function for some additional cleanup on odor descriptors
def clean_odors(odors):
    temp = odors.split(',')
    return ','.join(odor.strip().lower() for odor in temp)

# create dataframe for behavior
behavior = dfRaw.copy().drop(['Molecule Name'], axis=1).set_index('CID').sort_index()
behavior['Filtered Descriptors'] = behavior.apply(lambda row: clean_odors(row['Refined Descriptors']), axis=1)
behavior.drop(['Refined Descriptors'], axis=1).head()

# write to disk
molecules.to_csv('molecules.csv')
behavior.drop(['Refined Descriptors'], axis=1).to_csv('behavior.csv')