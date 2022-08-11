#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume
from pyrfume.odorants import get_cids, from_cids

cas = pd.read_csv('cas.txt', header=None)
cas = cas[0].apply(lambda x: x if '-' in x else None).dropna().values

cas_cid = get_cids(cas)

cas_cid['53896-26-7'] = 8892  # Hexanoic acid

cids = list(cas_cid.values())
molecules = pd.DataFrame(from_cids(cids)).set_index('CID').sort_index()

molecules = molecules[~molecules.index.duplicated()] # Remove duplicates
molecules.head()

weights_ = pd.read_csv('dragon-weights.csv')  # Extracted directly from supplement
symbols = pd.read_csv('dragon6-symbols.csv')  # Looked up in Dragon 6, for comparison to other versions

weights = weights_.join(symbols[['Dragon 6.0 symbol']]).set_index('Dragon 6.0 symbol')
weights = weights.drop('Descriptor', axis=1)
weights.head()

features = pd.read_csv('features.csv').set_index('PubChemID')

assert list(features.columns) == list(weights.index)

weighted_features = features.mul(weights['Weight'])
weighted_features.head()

weighted_features = features.T.mul(weights, axis=0)
weighted_features.head()

weighted_features[~weighted_features.isnull().all(axis=1)]

# Use API to compute distances

#from scipy.spatial.distance import pdist, squareform
#cids = features.index
#distances = pd.DataFrame(index=cids, columns=cids)
#distances[:] = squareform(pdist(weighted_features, metric='euclidean'))
#distances.to_csv('distances.csv')

# Write to disk
molecules.to_csv('molecules.csv')
weights.to_csv('weights.csv')
weighted_features.to_csv('features_weighted.csv')

