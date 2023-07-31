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

'''Benchmarking workflow for Arshamian 2022.

- Runs regression suite on Pleasantness (from behavior_1.csv) and Intensity (from behavior_2.csv) ratings
- Ratings are averaged across study participants (for Pleasantness) and molecule (for Intensity)
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
archive = 'arshamian_2022'
targets = ['Pleasantness', 'Intensity']
feature_sets = ['mordred', 'morgan']

# Pipelines to use
pipelines = [
    pbm.Model(estimator)
    for estimator in pbm.list_default_estimators('regression')
]


# -

# Function to prepare dataset; archive, prediction target, and feature set specific
def prepare_dataset(archive, target, feature_set):
    # Pleasantness ratings data
    if target == 'Pleasantness':
        # Stimulus ID = CID
        df = pyrfume.load_data(f'{archive}/behavior_1.csv').rename_axis('CID').reset_index()
        df.rename(columns={'Ranking': 'Pleasantness'}, inplace=True)
        # Average over participants
        df = df.groupby(by='CID').mean().drop(columns=['ParticipantID'])

    elif target == 'Intensity':
        df = pyrfume.load_data(f'{archive}/behavior_2.csv')
        df = pd.melt(df, var_name='Intensity', value_name='CID', ignore_index=False)
        df['Intensity'] = df['Intensity'].str.split(' ').str[-1].astype(int)
        # Average over CID
        df = df.groupby(by='CID').mean()
    
    # Convert to PyrfumeDataset class
    dataset = pbm.PyrfumeDataset(
        archive=archive,
        df=df,
        feature_set=feature_set,
        task='regression'
    )
    
    dataset.set_n_splits(5)
    
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
