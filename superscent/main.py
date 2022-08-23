#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pyrfume

# Load list of SMILES
cids = pd.read_csv('superscent_cids.txt', sep=' ', header=0)
cids.head()

# Get molecule info
molecules = pd.DataFrame(pyrfume.from_cids(cids.CID.to_list())).set_index('CID').sort_index()

molecules.head()

# Write to disk
molecules.to_csv('molecules.csv')

