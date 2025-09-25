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

'''Benchmarking workflow for aromadb.

- Runs classification tasks on 'Filtered Descriptors' from behavior.csv.
- Mordred and Morgan features sets are tried (independently, not merged).
- Descriptors appearing less than 5 times are dropped from the dataset.
- Train/test splits are generated using StratifiedKfold with n_splits=5.
- Train/test splits can be reporduced using indices returned by pyrfume.benchmarking.get_train_test_splits(dataset)
    where 'dataset' has been prepared using the prepare_dataset() function in this script.
'''

import pyrfume
import pandas as pd
import numpy as np
import pyrfume.benchmarking as pbm

# Dataset attributes
archive = 'aromadb'
targets = ['Filtered Descriptors'] 
feature_sets = ['mordred', 'morgan']

# Pipelines 
pipelines = [
    pbm.Model(estimator)
    for estimator in pbm.list_default_estimators('classification')
]


# Prepare the dataset
def prepare_dataset(archive, target, feature_set):
    # Load behavior data
    behavior = pyrfume.load_data(f'{archive}/behavior.csv').sort_index()
    behavior = behavior[~behavior.index.duplicated()]['Filtered Descriptors'].to_frame()

    # Retain only the first descriptor in multi-descriptor labels (these are 89/869 molecules)
    behavior['Filtered Descriptors'] = behavior['Filtered Descriptors'].str.split(';').str[0]
    
    # Stimulus column is CIDs
    df = behavior['Filtered Descriptors'].rename_axis('CID').to_frame()

    dataset = pbm.PyrfumeDataset(
        archive=archive,
        df=df,
        feature_set=feature_set,
        task='classification'
    )
    # Drop low frequency labels
    dataset.threshold_labels(min_counts=5)
    
    # Encode descriptors
    dataset.encode_labels()
    
    return dataset


pbm.verify_batch_settings(archive, targets=targets, feature_sets=feature_sets, prepare_dataset=prepare_dataset)

# +
kwargs = {'verbose' : 1}

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
