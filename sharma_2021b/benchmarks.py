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

'''Benchmarking workflow for sharma_2021b.

- Runs classification tasks on 'Primary Odor' and 'Sub Odor' from behavior.csv.
- Mordred and Morgan features sets are tried (independently, not merged).
- Descriptors appearing less than 10 times are dropped from the dataset.
- Train/test splits are generated using StratifiedKfold with n_splits=5.
- Train/test splits can be reproduced using indices returned by pyrfume.benchmarking.get_train_test_splits(dataset)
    where 'dataset' has been prepared using the prepare_dataset() function in this script.
'''

import pyrfume
import pandas as pd
from pyrfume import benchmarking as pbm

# +
# Dataset parameters
archive = 'sharma_2021b'
targets = ['Primary Odor', 'Sub Odor']
feature_sets = ['mordred', 'morgan']

# Pipelines to use
pipelines = [
    pbm.Model(estimator)
    for estimator in pbm.list_default_estimators('classification')]


# -

# Package dataset
def prepare_dataset(archive, target, feature_set):
    behavior = pyrfume.load_data(f'{archive}/behavior_1.csv').sort_index()
    behavior = behavior[~behavior.index.duplicated()] # Remove any duplicate molecules
    
    # Index is CID
    df = behavior[[target]].rename_axis('CID')

    # Convert to PyrfumeDataset class
    dataset = pbm.PyrfumeDataset(
        archive=archive,
        df=df,
        feature_set=feature_set,
        task='classification'
    )
    
    # Drop low frequency labels
    dataset.threshold_labels(min_counts=10)
    
    # Encode descriptors
    dataset.encode_labels()
    
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

# Dispaly the results
results.head()

# Display as heatmap
pbm.plot_heatmap(results)

# Compare to results from Dummy estimator
dataset = prepare_dataset(archive, targets[0], feature_sets[0])
dummy = pbm.evaluate_dummy_model(dataset)
dummy.head()

# Filter for top scores
best_results = pbm.get_best_results(results.copy(), include_pipeline_steps=True)
print(best_results.shape)
best_results.head(12)


# Save benchmarks
pbm.save_benchmarks(results, 'benchmarking.csv')

# Pickle best models, in addition to the 'prepare_dataset' function 
pbm.pickle_best_models(results=results, archive=archive, prepare_dataset=prepare_dataset)

# Execute the remote notebook for visualization
pbm.execute_viz_notebook(archive)
