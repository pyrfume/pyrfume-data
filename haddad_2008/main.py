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
from pyrfume.odorants import get_cids, from_cids

cas = pd.read_csv('cas.txt', header=None)
cas = cas[0].apply(lambda x: x if '-' in x else None).dropna().values

cas_cid = get_cids(cas)

cas_cid['53896-26-7'] = 8892  # Hexanoic acid

cids = list(cas_cid.values())
molecules = pd.DataFrame(from_cids(cids)).set_index('CID')
molecules.to_csv('molecules.csv')

weights_ = pd.read_csv('dragon-weights.csv')  # Extracted directly from supplement
symbols = pd.read_csv('dragon6-symbols.csv')  # Looked up in Dragon 6, for comparison to other versions

weights = weights_.join(symbols[['Dragon 6.0 symbol']]).set_index('Dragon 6.0 symbol')
weights = weights.drop('Descriptor', axis=1)
weights.to_csv('weights.csv')

features = pd.read_csv('features.csv').set_index('PubChemID')

assert list(features.columns) == list(weights.index)

weighted_features = features.mul(weights['Weight'])

weighted_features = features.T.mul(weights, axis=0)
weighted_features.to_csv('features_weighted.csv')

# +
# Use API to compute distances

#from scipy.spatial.distance import pdist, squareform
#cids = features.index
#distances = pd.DataFrame(index=cids, columns=cids)
#distances[:] = squareform(pdist(weighted_features, metric='euclidean'))
#distances.to_csv('distances.csv')
