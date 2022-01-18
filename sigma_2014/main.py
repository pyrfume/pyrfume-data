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

import numpy as np
import pandas as pd
import platform
import pyrfume
from pyrfume.odorants import from_cids, get_cids


# +
def get_data():
    with open("sigma_ff_catalog.txt", "r") as f:
        text = f.read()
        lines = text.split("\n")

    data = {}
    organoleptic = 0
    for line_num, line in enumerate(lines):
        if len(line):
            if not organoleptic and line[0] == "[":
                key = line.split("]")[0][1:]
                if platform.python_version() > "3.0":
                    key = key.replace("\u2011", "-")
                else:
                    key = key.decode("utf-8").replace("\u2011", "-").encode("ascii")
                organoleptic = 1
            if organoleptic and "Organoleptic" in line:
                try:
                    value = line.split(":")[1][1:]
                    if value[-1] in ["-", ","]:
                        if value[-1] == "-":
                            value = value[:-1]
                        else:
                            value = value + " "
                        value += lines[line_num + 1]
                    value = [i.strip() for i in value.split(";") if len(i.strip())]
                    data[key] = value
                    organoleptic = 0
                except Exception:
                    pass

    descriptors = []
    for x in data.values():
        descriptors += x
    descriptors = list(set(descriptors))  # Remove duplicates.
    return descriptors, data

descriptors, data = get_data()
data = pd.Series(data)

# + jupyter={"outputs_hidden": true} tags=[]
cas = data.index.tolist()
cids = get_cids(cas)
# -

data = data.to_frame('descriptors')
data.index.name = 'CAS'
data['CID'] = data.index.map(cids.get).astype(int)
data = data.reset_index().set_index('CID')

molecules = pd.DataFrame(from_cids(data.index)).set_index('CID').sort_index()
molecules.to_csv('molecules.csv')
molecules.head()

# Replace odorants with no CIDs with negative placeholders
ids = data[['CAS']].sort_index()
n_missing = (ids.index == 0).sum()
ids.index = np.arange(-1, -n_missing-1, -1).tolist() + ids.index[n_missing:].tolist()
ids = ids.sort_index()
ids.to_csv('identifiers.csv')
ids.head()

behavior = data.copy()
ids_reverse = ids.reset_index().set_index('CAS')
behavior.index = behavior['CAS'].apply(ids_reverse['index'].get)
behavior.index.name = 'CID'
behavior = behavior.drop('CAS', axis=1).sort_index()
behavior.to_csv('behavior_sparse.csv')

for descriptor in descriptors:
    behavior[descriptor] = behavior['descriptors'].apply(lambda x: descriptor in x)
behavior.drop('descriptors', axis=1).T.sort_index().T.astype(int)
behavior.to_csv('behavior.csv')
