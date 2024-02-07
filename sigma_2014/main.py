#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import platform
import pyrfume
from pyrfume.odorants import from_cids, get_cids

def get_data():
    with open("sigma_ff_catalog.txt", "r", encoding='utf-8') as f:
        text = f.read()
        lines = text.split("\n")

    data = {}
    count = 0
    for line_num, line in enumerate(lines):
        if not len(line):
            continue
        # getting the last CAS number
        if line.startswith("["):
            key = line.split("]")[0][1:]
            key = key.replace("\u2011", "-")

        # getting the descriptors
        if "Organoleptic:" in line:
            count += 1
            values = line.split(":")[1].strip()
            # if the line ends with a semicolon, it means the descriptor continues on the next line
            # if the line ends with a comma, it means the onion descriptor continues on the next line
            if values.endswith(";") or values.endswith(",") or values.endswith("wine-"):
                values += lines[line_num + 1].strip()
            # if the line ends with a dash, it means the descriptor has been split across two lines
            if values.endswith("-"):
                values = values[:-1] + lines[line_num + 1].strip()

            values = values.replace(" ", "", )
            valuelist = values.split(";")
            if key in data.keys():
                data[key].extend(valuelist)
                count -= 1
            else:
                data[key] = list(valuelist)

    #print(len(data))
    #print(count)
    #print(data['100-06-1'])
    if len(data) != count:
        raise ValueError("Not all molecules collected.")
    descriptors = set()
    for x in data.values():
        descriptors.update(x)
    descriptors = list(descriptors)
    return descriptors, data

descriptors, data = get_data()
data = pd.Series(data)

cas = data.index.tolist()
cids = get_cids(cas)

data = data.to_frame('descriptors')
data.index.name = 'CAS'
data['CID'] = data.index.map(cids.get).astype(int)
data = data.reset_index().set_index('CID')

molecules = pd.DataFrame(from_cids(data.index)).set_index('CID').sort_index()
molecules = molecules[~molecules.index.duplicated()]
molecules.to_csv('molecules.csv')
molecules.head()

# Replace odorants with no CIDs with negative placeholders
ids = data[['CAS']].sort_index()
n_missing = (ids.index == 0).sum()
ids.index = np.arange(-1, -n_missing-1, -1).tolist() + ids.index[n_missing:].tolist()
ids = ids.sort_index()
ids.index.name = 'Stimulus'
ids.to_csv('stimuli.csv')
ids.head()

behavior = data.copy()
ids_reverse = ids.reset_index().set_index('CAS')
behavior.index = behavior['CAS'].apply(ids_reverse['Stimulus'].get)
behavior.index.name = 'CID'
behavior = behavior.drop('CAS', axis=1).sort_index()
behavior.index.name = 'Stimulus'
behavior.to_csv('behavior_sparse.csv')
behavior.head()

for descriptor in descriptors:
    behavior[descriptor] = behavior['descriptors'].apply(lambda x: descriptor in x)
behavior.drop('descriptors', axis=1).T.sort_index().T.astype(int)
behavior.to_csv('behavior.csv')
behavior.head()
