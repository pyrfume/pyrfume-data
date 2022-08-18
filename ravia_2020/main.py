#!/usr/bin/env python
# coding: utf-8

# Ravia 2020 preprocessing
#
# Explanation of data
# Experiment 1: Similarity rating between several mixtures varying in intensity
#
# Experiment 2: Similarity rating between several mixtures varying in intensity
#
# Experiment 3 - Similarity ratings not included in Supplementary:
#   Bell challenge: similarity rating between three mixtures by perfumer
#
# Experiment 4: 50 pairs of mixture discriminated with triangle test by 27 subjects in replicate. 100 trials/subject
#
# Experiment 5: 50 pairs of mixtures discriminated with two-alternative same–different task by 30 subjects in replicate.
#   Each odorant was used at two different concentration. 400 trials/subject
#
# Experiment 6: 40 pairs of mixtures discriminated with two-alternative same–different task by 25 subjects in replicate.
#   Each odorant was used at two different concentration. 320 trials/subject
#
# Experiment 7 - Excluded from this analysis: Metemors discrimination on individual level

import pandas as pd
import numpy as np
import pyrfume

s1 = pd.read_csv('tableS1.csv')
s1.head()

# Compile all CIDs into list
mixture_columns = [col for col in s1.columns if 'Mix Component' in col]

all_cids = s1.CID.replace('Blank', np.nan).dropna().to_list()    + pd.melt(s1[mixture_columns].copy()).value.replace(0, np.nan).dropna().astype(int).to_list()

all_cids = list(set(all_cids))

# Get molecule info
molecules = pd.DataFrame(pyrfume.from_cids(all_cids)).set_index('CID').sort_index()

molecules.head()

# Create dataframe for stimuli.csv
stimuli = s1.drop([col for col in s1.columns if 'Unnamed:' in col], axis=1).copy()
stimuli[mixture_columns] = stimuli[mixture_columns].replace(0, np.nan)

# Merge mixture CIDs into CID columns
mixture_mask = (stimuli.Type != 'mono-molecule')
stimuli.loc[mixture_mask, 'CID'] = stimuli.loc[mixture_mask, mixture_columns]    .apply(lambda s: ';'.join(str(int(e)) for e in s if not pd.isna(e)), axis=1)

# Can now drop Mix Component columns and other columsn not relevant to stimuli description
drop_columns = ['Intensity1', 'Intensity2','Intensity3', 'PredictedDistance']
stimuli.drop(mixture_columns + drop_columns, axis=1, inplace=True)

# Each row is unique, so use row index as stimulus ID
stimuli.index.name = 'Stimulus'

stimuli.head()

# Intensity data -> behavior_1.csv
intensity = s1[(s1.Type == 'mono-molecule') & (s1.Experiment != 'Exp7')][['Intensity1', 'Intensity2', 'Intensity3']].copy()
intensity.index.name = 'Stimulus'

# Convert to long format
intensity = pd.melt(intensity, var_name='Dilution #', value_name='Intensity', ignore_index=False)
intensity['Dilution #'] = intensity['Dilution #'].str[-1]
intensity = intensity.reset_index().set_index(['Stimulus', 'Dilution #']).sort_index()
intensity.head()

s2 = pd.read_csv('tableS2.csv')
s2.head()

# Similarity data -> behavior_2.csv
similarity = s2.drop('Number', axis=1).copy()

# Convert mixture ID's to stimulus ID
mixtures1 = stimuli[(stimuli.Type != 'mono-molecule') & (stimuli.Experiment == 'Exp1')][['Experiment', 'ID']].copy()
mixtures1.ID = mixtures1.ID.str[-2:].astype(int)
mixtures2 = stimuli[(stimuli.Type != 'mono-molecule') & (stimuli.Experiment == 'Exp2')][['Experiment', 'ID']].copy()
mixtures2.ID = mixtures2.ID.str[-2:].astype(int)

mix_to_stim1 = dict(zip(mixtures1.ID, mixtures1.index))
mix_to_stim2 = dict(zip(mixtures2.ID, mixtures2.index))

similarity.loc[similarity.Experiment == 'Exp1', 'Stimulus 1'] = similarity.loc[similarity.Experiment == 'Exp1', 'Mixture1Id']     .map(mix_to_stim1)
similarity.loc[similarity.Experiment == 'Exp1', 'Stimulus 2'] = similarity.loc[similarity.Experiment == 'Exp1', 'Mixture2Id']     .map(mix_to_stim1)
similarity.loc[similarity.Experiment == 'Exp2', 'Stimulus 1'] = similarity.loc[similarity.Experiment == 'Exp2', 'Mixture1Id']     .map(mix_to_stim2)
similarity.loc[similarity.Experiment == 'Exp2', 'Stimulus 2'] = similarity.loc[similarity.Experiment == 'Exp2', 'Mixture2Id']     .map(mix_to_stim2)

similarity[['Stimulus 1', 'Stimulus 2']] = similarity[['Stimulus 1', 'Stimulus 2']].astype(int)
similarity = similarity.drop(['Experiment', 'Mixture1Id', 'Mixture2Id'], axis=1)     .set_index(['Stimulus 1', 'Stimulus 2']).sort_index()

similarity.head()

s3 = pd.read_csv('tableS3.csv')
s3.head()

# Discrimination data -> behavior_3.csv
discrim = s3[s3.Experiment != 'Exp7'].drop(['Mixture3', 'Mixture4', 'CompareIndex'], axis=1).copy()

# Convert Mixture # to Stimulus ID
mask = (stimuli.Type != 'mono-molecule') & ((stimuli.Experiment == 'Exp4,5') | (stimuli.Experiment == 'Exp6'))
id_to_stim = dict(zip(stimuli.loc[mask, 'ID'].str[1:].astype(int), stimuli.loc[mask].index))

discrim['Stimulus 1'] = discrim.Mixture1.map(id_to_stim)
discrim['Stimulus 2'] = discrim.Mixture2.map(id_to_stim)
discrim = discrim.drop(['Mixture1', 'Mixture2', 'Experiment'], axis=1)    .set_index(['Stimulus 1', 'Stimulus 2', 'Subject', 'TrialOrderInSession']).sort_index()
discrim.head()

# Write to disk
molecules.to_csv('molecules.csv')
stimuli.to_csv('stimuli.csv')
intensity.to_csv('behavior_1.csv')
similarity.to_csv('behavior_2.csv')
discrim.to_csv('behavior_3.csv')

