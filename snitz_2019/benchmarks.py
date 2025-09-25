# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

'''Benchmarking workflow for Snitz 2019.

- Runs regression suite on descpritor applicabiltiy, averaged over Subject (= study participant) responses in behavior_1.csv
- Mixtures are dropped; regressions run for only single molecules
- Mordred and Morgan features sets are tried (independently, not merged).
- Train/test splits are generated using Kfold with n_splits=5.
- Train/test splits can be reporduced using indices returned by pyrfume.benchmarking.get_train_test_splits(dataset)
    where 'dataset' has been prepared using the prepare_dataset() function in this script.
'''

import pyrfume
import pandas as pd
import numpy as np
import pyrfume.benchmarking as pbm

# +
# Dataset parameters
archive = 'snitz_2019'
targets = [
    'Pleasant', 'Volatile', 'WeakStrong', 'Bitter', 'Chemical', 'Spice', 'Sweet', 'FreshRotten',
    'CleanDirty', 'Burnt', 'Fruity', 'Develops', 'Garlic'
]
feature_sets = ['mordred', 'morgan']

# Pipelines to use
pipelines = [
    pbm.Model(estimator)
    for estimator in pbm.list_default_estimators('regression')
]


# -

# Function to prepare dataset; archive, prediction target, and feature set specific
def prepare_dataset(archive, target, feature_set):
    # Load behavior data; stimulus ID = CID
    behavior = pyrfume.load_data(f'{archive}/behavior_1.csv')
    
    # Get CID from stimuli.csv
    stimuli = pyrfume.load_data(f'{archive}/stimuli.csv')
    # Drop stimulus 'KO' (Musk) since listed (in manuscript Table 1) as same CAS# and CID as 'HG' (Galaxolide)
    stimuli.drop('KO', axis=0, inplace=True)
    
    # Add CID to behavior data
    df = behavior.join(stimuli)
    df.dropna(subset='CID', inplace=True) # Drops mixtures
    df.CID = df.CID.astype(int)
    
    # Average over Subject (= study participant)
    df = df.groupby(by='CID').mean()
    
    # Pick out particular prediction target
    df = df[[target]]
    
    # Convert to PyrfumeDataset class
    dataset = pbm.PyrfumeDataset(
        archive=archive,
        df=df,
        feature_set=feature_set,
        task='regression'
    )
    
    dataset.set_n_splits(10)
    
    return dataset


# +
# pbm.verify_batch_settings(archive, targets=targets, feature_sets=feature_sets, prepare_dataset=prepare_dataset)

# +
# Batch over feature set, prediciton targets, and methods list
kwargs = {'verbose': 1}

results = pbm.batch_gridsearchcv(
    archive=archive,
    targets=targets,
    feature_sets=feature_sets,
    pipelines=pipelines,
    prepare_dataset=prepare_dataset,
    **kwargs
)
# -

# Batch results
print(results.shape)
results.head()

# Filter for top scores
best_results = pbm.get_best_results(results)
print(best_results.shape)
best_results.head()

# Heatmap of results
pbm.plot_heatmap(best_results)

# Save benchmarks
pbm.save_benchmarks(results, 'benchmarks.csv')

# Pickle best models, in addition to the 'prepare_dataset' function 
pbm.pickle_best_models(results=results, archive=archive, prepare_dataset=prepare_dataset)

# Execute the remote notebook for visualization
pbm.execute_viz_notebook(archive)
