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

'''Benchmarking workflow for Mayhew 2022.

- Runs classification suite on 'odor.class' (odor vs odorless) from behavior_1.csv.
- Mordred and Morgan features sets are tried (independently, not merged).
- Train/test splits are generated using StratifiedKfold with n_splits=10.
- Train/test splits can be reporduced using indices returned by pyrfume.benchmarking.get_train_test_splits(dataset)
    where 'dataset' has been prepared using the prepare_dataset() function in this script.
'''

import pyrfume
import pandas as pd
import numpy as np
import pyrfume.benchmarking as pbm

# +
# Dataset parameters
archive = 'mayhew_2022'
targets = ['odor.class'] # Descriptor 1 is 'primary descriptor' acording to IFRA raw data
feature_sets = ['mordred', 'morgan']

# Pipelines to use
pipelines = [
    pbm.Model(estimator)
    for estimator in pbm.list_default_estimators('classification')
]


# -

# Function to prepare dataset; archive, prediction target, and feature set specific
def prepare_dataset(archive, target, feature_set):
    # Behavior data
    behavior = pyrfume.load_data(f'{archive}/behavior_1.csv')[[target]]
    # Get CIDs from stimuli
    stimuli = pyrfume.load_data(f'{archive}/stimuli.csv')
    # Add CID to prediction target
    df = behavior.join(stimuli).set_index('CID').sort_index()
    df = df[~df.index.duplicated()]

    # Convert to PyrfumeDataset class
    dataset = pbm.PyrfumeDataset(
        archive=archive,
        df=df,
        feature_set=feature_set,
        task='classification'
    )

    # Encode descriptors
    dataset.encode_labels()

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
