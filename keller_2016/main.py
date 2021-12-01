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

molecules.set_index('CID').to_csv('molecules.csv')

behavior.index = behavior.index.reorder_levels([1, 0, 2, 3]) # Put CID first

behavior = behavior.sort_index(level=0) # Sort by CID ascending

behavior.to_csv('behavior.csv')
