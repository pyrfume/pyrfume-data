# ---
# jupyter:
#   jupytext:
#     formats: py:light,ipynb
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

'''Benchmarking workflow for keller_2012.

- Runs regression tasks on 'intensity' and 'pleasantness' from behavior3.csv and behavior4.csv, respectively.
- Mordred and Morgan features sets are tried (independently, not merged).
- Train/test splits are generated using StratifiedKfold with n_splits=5.
- Train/test splits can be reproduced using indices returned by pyrfume.benchmarking.get_train_test_splits(dataset)
    where 'dataset' has been prepared using the prepare_dataset() function in this script.
'''

import pyrfume
from pyrfume import benchmarking as pbm
import pandas as pd
import numpy as np

# +
# Dataset parameters
archive = 'keller_2012'
targets = ['intensity', 'pleasantness'] # Not including 'threshold' as a target, because it only has 3 molecules
feature_sets = ['mordred', 'morgan']

# Pipelines to use
pipelines = [
    pbm.Model(estimator)
    for estimator in pbm.list_default_estimators('regression')
]

# Uncomment to use this reduced pipeline for testing purposes
# estimators = ['Ridge', 'Lasso']
# reduced_param_grids = [
#     {"alpha": [0.01, 0.1]},
#     {"alpha": [0.01, 0.1], "max_iter": [1000, 3000, 5000, 7000]}
# ]
# pipelines = [
#     pbm.Model(estimator, param_grid=param_grid)
#     for estimator, param_grid in zip(estimators, reduced_param_grids)
# ]
# -

# Function to prepare dataset; archive, prediction target, and feature set specific
def prepare_dataset(archive, target, feature_set):
    # Subset of stimuli for the intensity (behavior_3) and pleasantness (behavior_4) data
    stimuli = pyrfume.load_data(f'{archive}/stimuli.csv')
    stimuli['experiment type'] = stimuli.index.str[0]
    
    if target == 'intensity':
        behavior = pyrfume.load_data(f'{archive}/behavior_3.csv')
        stimuli = stimuli[(stimuli['experiment type'] == 'I') & (stimuli['concentration'] == 'high')]['CID'].to_frame()
        stimuli = stimuli.dropna(axis=0) 
   
    elif target == 'pleasantness':
        behavior = pyrfume.load_data(f'{archive}/behavior_4.csv')
        stimuli = stimuli[(stimuli['experiment type'] == 'I') & (stimuli['concentration'] == 'high')]['CID'].to_frame()
        stimuli = stimuli.dropna(axis=0)
    else:
        raise ValueError(f'Invalid target: {target}')
            
    # Aggregate by stimulus, calculate mean target value of target variable, and merge in CIDs
    mean_target = behavior.groupby(by='stimulus').mean()[target].to_frame()
    df = pd.merge(stimuli, mean_target, left_index=True, right_index=True).set_index('CID').sort_index()
    df.index = df.index.astype(int)
    
    #Convert to PyrfumeDataset class
    dataset = pbm.PyrfumeDataset(
        archive=archive,
        df=df,
        feature_set=feature_set,
        task='regression'
    )
    
    return dataset


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
best_results = pbm.get_best_results(results.copy(), include_pipeline_steps=True)
print(best_results.shape)
best_results.head(12)

# Heatmap of results
pbm.plot_heatmap(results)

# Compare best score results to those from Dummy estimators
dataset = prepare_dataset(archive, targets[0], feature_sets[0])
dummy = pbm.evaluate_dummy_model(dataset)
dummy.head()

# Save benchmarks
pbm.save_benchmarks(results, 'benchmarking.csv')

# +
# Pickle best models, in addition to the 'prepare_dataset' function 
pbm.pickle_best_models(results=results, archive=archive, prepare_dataset=prepare_dataset)

#Execute the remote notebook for visualization
pbm.execute_viz_notebook(archive)
