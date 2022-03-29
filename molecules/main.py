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

from itertools import chain
import pandas as pd
import pyrfume
from pyrfume.odorants import from_cids

archives = pyrfume.list_archives()
archive_cids = {}
for archive in archives:
    if archive not in ['molecules']:
        molecules = pyrfume.load_data(f'{archive}/molecules.csv')
        archive_cids[archive] = molecules.index
cids = sorted(chain(*archive_cids.values()))

usage = pd.DataFrame(index=cids)
for archive in archives:
    usage[archive] = usage.index.isin(archive_cids[archive]).astype(int)
usage.to_csv('usage.csv')

molecules = pd.DataFrame(from_cids(cids)).set_index('CID')
molecules.to_csv('molecules.csv')


