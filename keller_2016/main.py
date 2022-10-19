#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
from pyrfume import keller
from pyrfume.odorants import from_cids

# Loads the raw data from 12868_2016_287_MOESM1_ESM.xlsx
raw = keller.load_raw_bmc_data()
raw.head()

behavior = keller.format_bmc_data(raw,
                                  only_dream_subjects=False,  # Whether to only keep DREAM subjects
                                  only_dream_descriptors=False,  # Whether to only keep DREAM descriptors
                                  only_dream_molecules=False)  # Whether to only keep DREAM molecules))

behavior.head()

cids = behavior.index.unique(level='CID')
molecules = pd.DataFrame(from_cids(cids))
molecules = molecules.set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()] # Remove duplicates if any
molecules.head()

behavior.index = behavior.index.reorder_levels([1, 0, 2, 3]) # Put CID first
behavior = behavior.sort_index(level=0) # Sort by CID ascending
behavior.index.rename('Stimulus', level='CID', inplace=True)
behavior.head()

# Create dataframe for stimuli.csv; all stimului are CIDs
stimuli = pd.DataFrame(molecules.index.copy(), index=molecules.index.copy())
stimuli.index.name = 'Stimulus'
stimuli.head()

# Write to file
behavior.to_csv('behavior.csv')
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')

