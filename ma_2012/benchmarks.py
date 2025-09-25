# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

'''Benchmarking workflow for ma_2012.

- Runs classification tasks on 'Dorsal Response from behavior1.csv.
- Mordred and Morgan features sets are tried (independently, not merged).
- Train/test splits are generated using StratifiedKfold with n_splits=5.
- Train/test splits can be reporduced using indices returned by pyrfume.benchmarking.get_train_test_splits(dataset)
    where 'dataset' has been prepared using the prepare_dataset() function in this script.
'''

import pyrfume
import pandas as pd
import numpy as np
import benchmarking as pbm

# +
# Dataset parameters
archive = 'ma_2012'
targets = ['Dorsal Response']
feature_sets = ['mordred', 'morgan'] # mordred, morgan, or mordred_morgan (to use both)

# Pipelines to use
pipelines = [
    pbm.Model(estimator)
    for estimator in pbm.list_default_estimators('classification')
]

# Code below is for a testing on a reduced set. Comment-out or delete to run the real pipeline
# estimators = ['DecisionTreeClassifier', 'RandomForestClassifier', 'BernoulliNB']
# reduced_param_grids = [
#     {'min_impurity_decrease': [0.001, 0.005], 'max_features': ['sqrt', 'log2']},
#     {'n_estimators': [100], 'min_impurity_decrease': [0.001, 0.005], 'max_features': ['sqrt', 'log2']},
#     {'alpha': [1, 10], 'fit_prior': [True], 'binarize': [0.25, 0.5, 0.75]}
# ]

# pipelines = [
#     pbm.Model(estimator, param_grid=param_grid)
#     for estimator, param_grid in zip(estimators, reduced_param_grids)
# ]
# -

# Function to prepare dataset; archive, prediction target, and feature set specific
# This part is idiosyncratic to each archive
def prepare_dataset(archive, target, feature_set):
    # Load behavior data
    df = pyrfume.load_data(f'{archive}/behavior_1.csv')

    # Stimulus ID is just CID
    df.index.name = 'CID'

    # Map DR integer values back to meaning presented by Ma et al
    df[target] = df[target].map({0: 'None', 1: 'Low', 2: 'Medium', 3: 'High'})
    
    # Convert to PyrfumeDataset class
    dataset = pbm.PyrfumeDataset(
        archive=archive,
        df=df,
        feature_set=feature_set,
        task='classification'
    )
    
    # Encode DR labels
    dataset.encode_labels()

    return dataset


prepare_dataset(archive, targets[0], feature_sets[0])

# +
# Batch process over feature set, prediciton targets, and methods list
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

# Batch results for best scores
print(results.shape)
results.head()

# Heatmap of results
pbm.plot_heatmap(results)

# Filter for top scores
best_results = pbm.get_best_results(results.copy(), include_pipeline_steps=True)
print(best_results.shape)
best_results.head(12)

# Save benchmarks
pbm.save_benchmarks(results, 'benchmarking.csv')

# Pickle best models, in addition to the 'prepare_dataset' function 
pbm.pickle_best_models(results=results, archive=archive, prepare_dataset=prepare_dataset)

# Execute the remote notebook for visualization
pbm.execute_viz_notebook('viz_template_notebook.ipynb', f'{archive}_viz.ipynb', archive)
