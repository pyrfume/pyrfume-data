#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import os
import requests
import bs4

# Scrape functions
def get_codebook(url):
    resp = requests.get(url) 
    
    soup = bs4.BeautifulSoup(resp.content)
    return soup.find('div', {'id': 'Codebook'})

def codebook_to_answer_key(codebook, question):
    if 'CSQ260' in question: # On website sub-questions for CSQ260 have lower case letters
        question = question[:-1] + question[-1].lower()
        
    answer_key = {}
    table = codebook.find('h3', {'id': question}).findNext('tbody')
    for tr in table.find_all('tr'):
        row = tr.text.strip().split('\n')
        if row[0] == '.':
            answer_key['-1'] = row[1]
        elif 'to' not in row[0]:
            answer_key[row[0]] = row[1]
    return answer_key

def codebook_to_meta_data(codebook, question):
    if 'CSQ260' in question: # On website sub-questions for CSQ260 have lower case letters
        question = question[:-1] + question[-1].lower()
        
    tmp = codebook.find('h3', {'id': question}).findNext('dl')
    sas_label = tmp.find('dt', text='SAS Label: ').findNext('dd').text.strip()
    english_text = tmp.find('dt', text='English Text: ').findNext('dd').text.strip()
    target = tmp.find('dt', text='Target: ')        .findNext('dd').text.strip().replace('\r','').replace('\n','').replace('\t','').replace(' -', ' - ')
    return [sas_label, english_text, target]

# "The Taste and Smell Questionnaire Section (variable name prefix CSQ) collected interview data on self-reported taste 
# and smell ability, selected symptoms of and medical treatment for taste and smell disorders, and data on conditions that
# may represent risk factors for taste and smell disorders. The CSQ questionnaire was designed to provide data to support
# the Healthy People 2020 objectives for taste and smell disorders (Healthy People, 2020)."
#
# https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/CSQ_H.htm

csq = pd.read_sas('CSQ_H.XPT').fillna(-1).astype(int).astype(str)
csq.rename(columns={'SEQN': 'Subject'}, inplace=True)
csq.Subject = csq.Subject.astype(int)

print(csq.shape)
csq.head()

# Scrape to get questionairre meta data and answer keys
csq_codebook = get_codebook('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/CSQ_H.htm')

csq_meta_data = {}
for question in [col for col in csq.columns if col != 'Subject']:
    csq_meta_data[question] = codebook_to_meta_data(csq_codebook, question)
    
    # Map numerical answers to text answers for each question
    answer_key = codebook_to_answer_key(csq_codebook, question)
    csq[question].replace(answer_key, inplace=True)
    
csq.head()

# Convert to long form for behavior_1.csv
behav1 = pd.melt(csq, id_vars='Subject', var_name='Stimulus', value_name='Response')
behav1 = behav1.set_index(['Stimulus', 'Subject']).sort_index()
behav1.head()

# "The taste and smell exam measured the ability to taste and smell, using an odor identification test and salt and quinine
# taste testing. The objectives of this component were:
#
# 1. to provide reference data for taste and smell testing for U.S. adults aged 40 and over;
# 2. to examine variations in the ability to smell and to taste salt and bitter tastants and analyze these variations with
# NHANES hypertension, nutritional, and obesity data; and
# 3. to help estimate the prevalence of U.S. adults who may not recognize the odor of smoke and natural gas, which are
# important early warning signals for home safety hazards."
#
# https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/CSX_H.htm
csx = pd.read_sas('CSX_H.XPT').fillna(-1)
csx.rename(columns={'SEQN': 'Subject'}, inplace=True)
csx.Subject = csx.Subject.astype(int)
csx['CSXTSEQ'] = csx['CSXTSEQ'].str.decode("utf-8")

for question in [col for col in csx.columns if col not in ['Subject', 'CSXTSEQ']]:
    csx[question] = csx[question].astype(int).astype(str)

print(csx.shape)
csx.head()

# Scrape to get questionairre meta data and answer keys
csx_codebook = get_codebook('https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/CSX_H.htm')

csx_meta_data = {}
for question in [col for col in csx.columns if col != 'Subject']:
    csx_meta_data[question] = codebook_to_meta_data(csx_codebook, question)
    
    # Map numerical answers to text answers for each question
    answer_key = codebook_to_answer_key(csx_codebook, question)
    csx[question].replace(answer_key, inplace=True)
    
csx.head()

# Convert to long form for behavior_2.csv
behav2 = pd.melt(csx, id_vars='Subject', var_name='Stimulus', value_name='Response')
behav2 = behav2.set_index(['Stimulus', 'Subject']).sort_index()
behav2.head()

# Combine meta data to create stimuli.csv; stimulus is question code
csq_meta = pd.DataFrame.from_dict(csq_meta_data, orient='index', columns=['SAS Label', 'English Text', 'Target'])
csq_meta.index.name = 'Stimulus'

csx_meta = pd.DataFrame.from_dict(csx_meta_data, orient='index', columns=['SAS Label', 'English Text', 'Target'])
csx_meta.index.name = 'Stimulus'

stimuli = pd.concat([csq_meta, csx_meta], axis=0)
stimuli.head()

# Write to disk
stimuli.to_csv('stimuli.csv')
behav1.to_csv('behavior_1.csv')
behav2.to_csv('behavior_2.csv')

