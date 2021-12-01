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

# # Ravia 2020 preprocessing
# ### See main.rmd for most of the pre-processing

import pandas as pd
from pyrfume.odorants import from_cids 
mixtures = pd.read_csv('mixtures.csv')

all_cids = set()
for cids in mixtures['Mixture.Cids']:
    cids = [int(x) for x in cids.split(' ')]
    all_cids |= set(cids)
all_cids = sorted(list(all_cids))

molecules = from_cids(all_cids)
molecules = pd.DataFrame(molecules).set_index('CID')

molecules.to_csv('molecules.csv')
