#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import numpy as np
import matplotlib
from skimage import data,io,viewer,feature,exposure
from pyrfume import odorants
from collections import OrderedDict


# ### First pass at grabbing molecules
#  - Grab all .png files from the directory of glomerular maps
#  - Remove suffixes from file names that are idiosyncratic to Leon
#  - Grab the PCIDs we can, from this first naive pass
# ***

# In[2]:


saved_path = os.getcwd()
root = '/Users/jcastro/Dropbox/CnG/' # set as needed
ccc = 'Leon/Glomerular archive/ccc copy' # dir of glomerular maps

os.chdir(root + ccc)
dir_list = [f for f in os.listdir() if f.endswith('.png')]
os.chdir(saved_path)

# format for pubchem search
filenames = [chem.split('.png')[0].replace('_',' ')              for chem in dir_list]

# drop irrelevant suffixes 
chemnames = [' '.join(chem.split()[:-1])              for chem in filenames] 

# filter out empty list entries
del_indices = [i for i in range(len(chemnames)) if chemnames[i] == '']
del_indices = sorted(del_indices, reverse=True)

for d in del_indices:
    chemnames.pop(d)
    filenames.pop(d)
    dir_list.pop(d)


# ***
# If you want to get a sense of the relationship btw chemnames, filenames. There are some pretty funky/idiosyncratic suffixes on the original file names:
# ***

# In[3]:


i = 420 # choose from: 0-539
print(chemnames[i], filenames[i], dir_list[i] )


# ***
# Grab PCIDs for all chemicals in the list. Use an ordered dict to maintain alignment w/ file names
# ***

# In[4]:


molecules = OrderedDict()
molecules = odorants.get_cids(chemnames)


# In[5]:


molecule_names = list(molecules.keys())


# ### Hunt down orphan molecules

# In[6]:


molecules_saved = molecules # just for safe-keeping
orphans = [key for key, value in molecules.items() if value == 0]

print("%2d molecules found, and %2d orphans" % (len(molecules)-len(orphans),len(orphans)))


# ***
# Dict of hand-searched (Google + Pubchem) orphans, including replacement names + PCIDs where possible
# ***

# In[7]:


# formatted as:  'orphan name' : ('replacement name', PCID)

orphan_lookup = {
   
    'pentanol fgrp' : ('pentanol',6276),
    'isobutyricacid' : ('isobutyric acid',6590),
    'acetone 45' : ('acetone',180),
    'a phe mint 1' : 0,
    'dietsuberate' : ('deithyl suberate',16301),
    '7 oxabicycloheptane' : ('Oxabicycloheptane',18675856),
    '2 nonanone kecon' : ('2-Nonanone',13187),
    'a ionone a io 1': 0,
    'trans 2 penten acid': 0,
    '5hydroxymethylfurfural' : ('5-Hydroxymethylfurfural',237332),
    '4 hydroxybenzal' : ('4-Hydroxybenzaldehyde',126),
    'isopropanol salcohols' : ('Isopropanol',3776),
    'methanol methcon': 0,
    'methyl3 methylbutenoa': 0,
    'a angelicalactone' : ('alpha-Angelica lactone',11559),
    'menicotinate' : ('Methyl nicotinate',7151),
    'dietmalonate' : ('Diethyl malonate',7761),
    'm anisaldehyde' : ('3-Methoxybenzaldehyde',11569),
    'menthone tc' : ('Menthone',26447),
    'pentadecane SG' : ('Pentadecane',12391),
    '2,3,5': 0,
    'menthylacatat mint 1' : ('Menthyl acetate',27867),
    'limonene( )' : ('Limonene',22311),
    'p tolylphenylacetate' : ('p-Tolyl phenylacetate',60997),
    '5methylfurfural' : ('5-Methylfurfural',12097),
    'pentanal fgrp' : ('Pentanal',8063),
    '5678tetrahydroquinoxaline' : ('5,6,7,8-Tetrahydroquinoxaline',36822),
    '3(5 methyl 2': 0,
    'terp4ol( )': 0,
    'ethyl butyrate etbuENR' : ('Ethyl butyrate',7762),
    'octanoicacid' : ('Octanoic acid',379),
    '25dime 4 hydr 3fura': 0,
    'ethylbenzene ebco' : ('Ethylbenzene',7500),
    'limonene(-)' : ('(-)-Limonene',439250),
    'methyloctanoate acros' : ('Methyl octanoate',8091),
    '4oxoisophorone' : ('2,6,6-Trimethyl-2-cyclohexene-1,4-dione',62374),
    'heptane 800' : ('Heptane',8900),
    'methylsalicylate MScon1' : ('Methyl salicylate',4133),
    'o anisaldehyde' : ('2-Methoxybenzaldehyde',8658),
    'benzylacohol' : ('Benzyl Alcohol',244),
    'methylval fgrp': 0,
    '2 hexanone fgrp' : ('2-Hexanone',11583),
    'methylthiosalicylate' : ('Methyl 2-mercapto-4-methylbenzoate',57555490),
    '4 pentenoicacid' : ('4-Pentenoic acid',61138),
    'methyl trans 2': 0,
    '2isobut3methoxpyrazine' : ('2-Isobutyle-3-methoxypyrazine',32594),
    'cyclohexanehigh' : ('Cyclohexane',8078),
    'mesalicylate MScon' : ('Methyl salicylate',4133),
    '2,3,5,6 tetramethylpyr' : ('2,3,5,6-Tetramethylpyridine',19562),
    'dipropyldisulfide' : ('Dipropyl disulfide',12377),
    'trans trans 26 nonadien': 0,
    '2 octanone kecon' : ('2-Octanone',8093),
    'cycloheptanelow' : ('Cycloheptane',9265),
    'phenylbenzoate' : ('Phenyl Penzoate',7169),
    'limonene(+)' : ('D-Limonene',440917),
    'valericacid vcn1' : ('Valeric acid',7991),
    '2 ketobutyricacid' : ('2-Oxobutanoic acid',58),
    'terp4ol(+)' : ('(+)-Terpinen-4-ol',2724161),
    'nerol cis' : ('Nerol',643820),
    'methyloctanoate aldrich' : ('Methyl octanoate',8091),
    '2 mecycpropcarbacid': 0,
    'ethyllactate' : ('Ethyl lactate',7344),
    '3acetyl25dimethylfuran':('3-Acetyl-2,5-dimethylfuran',61527),
    'caproicacid' : ('Hexanoic acid',8892),
    'methyl2furoate' : ('Methyl 2-furoate',11902),
    '1 pentanol fgrp' : ('1-Pentanol',6276),
    'cyclooctanehigh' : ('Cyclooctane',9266),
    'ethylcaproate' : ('Ethyl hexanoate',31265),
    'ethanol salcohols' : ('Ethanol',702),
    'menthylisoval mint 1' : ('1-menthyl isovalerate',129773535),
    'geraniol cis' : ('Geraniol',637566),
    'methylisocaproate' : ('Methyl 4-methylvalerate',17008),
    'lauricacid' : ('Lauric acid',3893),
    'alpha cedrene EO1' : ('alpha-Cedrene',442348),
    'methylcyclohexanehigh' : ('Methylcyclohexane',7962),
    'D(+) camphor' : ('D-Camphor',159055),
    'L menthone enr' : ('Menthone',26447),
    'decanal etbuENR' : ('Decanal',8175),
    'valeric acid fgrp' : ('Valeric acid',7991),
    'eucalyptol SG' : ('Eucalyptol',2758),
    'Leon': 0,
    'heptane 2500' : ('Heptane',8900),
    'cis jasmone' : ('Jasmone',1549018),
    'trans 2 tridecenal' : ('3-(Dimethylamino) acrylaldehyde',638320),
    'pentadecane alfaaesar' : ('Pentadecane',12391),
    'trans,trans 2,4 decadienal' : ('2,4-Decadienal',5283349),
    'geranylacetate cis' : ('Geranyl acetate',1549026),
    '4 tertbutylpyridine' : ('4-tert-Butylpyridine',19878),
    '2 butanone' : ('Methyl ethyl ketone',6569),
    'dietoxalate' : ('Diethyl oxalate',7268),
    'me anthranilate' : ('N-Methylanthranilate',6942392),
    'amyl acetate etbuENR': 0,
    'caproic acid aci1' : ('Hexanoic acid',8892),
    'menthone mint 1' : ('Menthone',26447),
    'propanol simp' : ('Propanol',1031),
    'cyclopentcarb acid' : ('Cyclopentanecarboxylic acid',18840),
    'pentadecane sigma' : ('Pentadecane',12391),
    'aromadendrene EO1' : ('Aromadendrene',91354),
    'butyricacid' : ('Butyric acid',264),
    'menthylacetat mint 1' : ('Menthyl acetate',27867),
    'heptane 250' : ('Heptane',8900),
    'trimecyclohexanone': 0,
    'caproicacid aci1' : ('Hexanoic acid',8892),
    'sorbicalsohol' : ('2,4-Hexadien-1-ol',641256),
    'methyl2octynoate' : ('Methyl 2-octynoate',8092),
    'L carvone' : ('(-)-Carvone',439570),
    '2,3,5 trimethylpyrazin' : ('2,3,5-Trimethylpyrazine',26808),
    '3,7 dimethyloctan 1': 0,
    'pentadecane acros' : ('Pentadecane',12391),
    'methyloctanoate' : ('Methyl octanoaate',8091),
    't butylacetic acid': 0,
    'delta dodecanalactone' : ('delta-Dodecalactone',12844),
    'undecylenicacid' : ('Undecylenic acid',5634),
    '2 butEtOHacetate' : ('2-Butoxyethyl acetate',8160),
    'me salicylate' : ('Methyl salicylate',4133),
    'menthol mint 1': 0,
    'me phenylacetate' : ('Methyl phenylacetate',7559),
    'pheanthranilate': 0,
    'cyclo octane' : ('Cyclooctane',9266),
    'aceticacid' : ('Acetic acid' ,176),
    'ethyltetradecanoate' : ('Ethyl tetradecanoate',31283),
    'cyclohexcarboxaldehyde' : ('Cyclohexanecarboxaldehyde',16275),
    'methl 2 octynoate' : ('Methyl 2-octynoate',8092),
    'trimethylaceticacid' : ('Pivalic acid',6417),
    'decanal DeTimCo' : ('Decanal',8175),
    '2acetylpyrazine' : ('Acetylpyrazine',30914),
    '2,6 dimethyl 5': 0,
    'phenethylhexanoate' : ('2-Phenylethyl hexanoate',61384),
    'dimethyloctenone' : ('4,7-Dimethyloct-6-en-3-one',102844),
    'limonene mint 1' : ('Limonene',22311),
    'celeriax': 0,
    'methylcyclohexanecarboxy': 0,
    'acetone tc' : ('Acetone',180),
    'methylo toluate' : ('Methyl o-toluate',33094),
    'ketonemoschus' : ('Musk ketone',6669),
    '14dimethoxybenzene' : ('1,4-Dimethoxybenzene',9016),
    'sec butanol' : ('2-Butanol',6568),
    'isoamylacetate' : ('Isoamyl acetate',31276),
    'trans 2 hexen 1 ol' : ('2-Hexen-1-ol',5318042),
    'methyltrans2octenoate' : ('Methyl oct-2-enoate',5364532),
    'methyl3 aminobenzoate' : ('Methyl 3-aminobenzoate',78274),
    'cyclamenaldehyde' : ('3-(4-Isopropylphenyl)-2-methylpropanal',517827),
    'trans 2 octene' : ('2-Octene',5364448),
    '2 decanone SG' : ('2-Decanone',12741),
    'trans2hexenoicacid' : ('trans-2-Hexanoic acid',5282707),
    'heptane 7500' : ('Heptane',8900),
    'hexylacetate' : ('Hexyl acetate',8908),
    'methylcyclohexanelow' : ('Methylcyclohexane',7962),
    'tertbutylbenzene' : ('Tert-butylbenzene',7366),
    '2 methylbutyric acid' : ('2-Methylbutanoic acid',8314),
    'p tolualdehyde' : ('4-Methylbenzaldehyde',7725),
    'cyclooctanelow' : ('Cyclooctane',9266),
    '4 methylvaleric acid' : ('4-Methylpentanoic acid',12587),
    'cinnamaldehyde SG' : ('Cinnamaldehyde',637511),
    'perrilylalcohol' : ('Perillyl alcohol',10819),
    'methylglycolate' : ('Methyl glycolate',66774),
    'monoethylsuccinate': 0,
    'methylvalerate fgrp' : ('Methyl valerate',12206),
    'mecinnamate' : ('Methyl cinnamate',637520),
    'methyl3 hydroxybenzoate' : ('Methyl 3-hydroxybenzoate',88068),
    'tiglicaldehyde' : ('Tiglic aldehyde',5321950),
    'trans 3 hexen acid': 0,
    'nerylacetate cis' : ('Neryl acetate',1549025),
    'a ionone a io enr': 0,
    '.ipynb': 0,
    'propylbutyrate' : ('Propyl butyrate',7770),
    'cyclohexanelow' : ('Cyclohexane',8078),
    'bornyl acetate etbuENR' : ('Bornyl acetate',6448),
    'valericacid fgrp' : ('Valeric acid',7991),
    'strawberryaldehyde' : ('Ethyl methylphenylglycidate',6501),
    'thujone SG24' : ('Thujone',261491),
    'valericacid enr' : ('Valeric acid',7991),
    'pentadecane fluka' : ('Pentadecane',12391),
    'o tolualdehyde' : ('2-Methylbenzaldehyde',10722),
    'cyclobutcarb acid': 0,
    'laurylacetate' : ('Dodecyl acetate',8205),
    '245trimethyloxazole' : ('2,4,5-Trimethyloxazole',30215),
    'carvone(-) enr' : ('(-)-Carvone',439570),
    'dietsuccinate': 0,
    '46dimethyl2pyrone' : ('4,6-Dimethyl-2H-pyran-2-one',12662),
    '4 methylnonanoicacid' : ('4-Methylnonanoic acid',62003),
    'allyl2furoate' : ('Allyl 2-furoate',61337),
    'trans cis 26 nonadien': 0,
    'methyloctanoate alfaesr' : ('Methyl octanoate',8091),
    '45dime3hydrox25dihyfur': 0,
    'dimetmalonate' : ('Dimethyl malonate',7943),
    'propanol salcohols' : ('Propanol',1031),
    'isobutylangelate' : ('Isobutyl angelate',5367807),
    '2octynoicacid' : ('2-Octynoic acid',21872),
    'acetone ac tc' : ('Acetone',180),
    'heptane 80' : ('Heptane',8900),
    'cycloheptanehigh' : ('Cycloheptane',9265),
    'isoamyltiglate' : ('Isoamyl tiglate',5463682),
    '2 me 3 buten 2 ol' : ('2-Methyl-3-buten-2-ol',8257),
    'methanol salcohols' : ('Methanol',887),
    '2,3,5,6': 0,
    'methylacetate' : ('Methyl acetate',6584),
    '5methylquinoxaline' : ('5-Methylquinoxaline',61670),
    'm tolualdehyde' : ('3-Methylbenzaldehyde',12105),
    'trans trans 26' : 0
} 


# ***
# It's useful to have a reversed lookup table of the orphans, so we can get original Leon chemical names (which failed to generate a PCID) from the new names:
# ***

# In[10]:


o_target = [] # orphan 'targets': new names of renamed molecules
o_source = [] # orphan 'sources': original names of molecules that were renamed

orphan_reassignments = list(orphan_lookup.values())
orphan_originals = list(orphan_lookup.keys())

for i in range(len(orphan_reassignments)):
    if orphan_reassignments[i] != 0:
        o_target.append(orphan_reassignments[i][0])
        o_source.append(orphan_originals[i])
        
orphans_reversed = OrderedDict()
orphans_reversed = {o_target[i] : o_source[i] for i in range(len(o_target))}


# ### Create the {new molecules : PCID} dictionary
# 
#  - dictionary containing the newly de-orphaned molecules

# In[11]:


# write a new 'molecules' dict, w renamed/de-orphaned molecules

molecules_new = OrderedDict() # dict to store the renamed molecules

for key, value in molecules.items():
    # does the molecule (in the original list) lack a PCID? 
    if value == 0: 
        # if the orphaned molecule was renamed & we then found its PCID:
        if (key in orphan_lookup) and orphan_lookup[key] != 0:
            new_name, newPCID = orphan_lookup[key]
            molecules_new[new_name] = newPCID # write the new PCID to the dict
        # otherwise, it remains orphaned
        else:
            molecules_new[key] = 0 # we couldn't find the PCID
    else:
        molecules_new[key] = value # original molecule was well-named 
        
molecule_names_new = list(molecules_new.keys())


# ### Write the file lookup table
# 
#  - Each row is a molecule name
#  - Each entry is a list (variable #s) of associated .png files 

# In[39]:


# store all files associated w/ a given molecule as an array

file_lookup = [] # final list of all .png images, aligned to the list of molecule names

cnames = np.array(chemnames)
filenames = np.array(dir_list)

for i in range(len(molecule_names_new)):
    idx_renamed = []
    idx_orig = []
    
    if molecule_names_new[i] in o_target:
        tname = orphans_reversed[molecule_names_new[i]]
        idx_renamed = np.where(cnames == tname)[0]
        
    idx_orig = np.where(cnames == molecule_names_new[i])[0]
    combined_indices = [i for i in idx_orig] + [i for i in idx_renamed]
    file_lookup.append(filenames[combined_indices])


# In[13]:


# A couple quick checks that # molecules & # rows in file lookup table match: 

print("There are %2d molecules in the molecule list, and %2d rows in the file lookup table \n" %       (len(molecule_names_new),len(file_lookup)))

# inspect individual molecules to check alignment btw. files and chem names

mol_number = 110 # (btw. 1-408)

print ('Molecule #: ',mol_number)
print('name: ', molecule_names_new[mol_number])
print('associated files: ',file_lookup[mol_number])


# In[14]:


# Quick and dirty summary of the data set:

images_per_molecule = []
max_images = 1

for i in range(len(file_lookup)):
    images_per_molecule.append(len(file_lookup[i]))

imagevals = list(range(1,max(images_per_molecule)+1))

for i in imagevals: 
    numfiles = len([n for n in images_per_molecule if n == i])
    print('There are %2d molecules with %2d associated files' % (numfiles, i))
    


# ### Write the csv files for pyrfume
# 

# In[44]:


PCIDs = [value for key,value in molecules_new.items()]

print('there are %2d PCIDs and %2d lines in the file lookup table'       % (len(PCIDs),len(file_lookup)))


# ***
# Remove all entries where PCIDs = 0 

# In[46]:


final_mols = PCIDs
final_files = file_lookup

PCID_zero_idx = sorted(np.where(np.array(PCIDs) == 0)[0], reverse=True)

for z in PCID_zero_idx:
    final_mols.pop(z) # final
    final_files.pop(z)

print(len(final_mols),len(final_files),len(file_lookup))


# In[47]:


mols = odorants.from_cids(PCIDs)
mols_df = pd.DataFrame(mols)
mols_df.set_index('CID', inplace = True)


# In[48]:


mols_df.head()


# In[49]:


mols_df.to_csv('molecules.csv')


# In[52]:


import csv

# a bit hack-y, but can't get each row to write properly otherwise:
csv_file = open('behavior.csv', 'w')
writer = csv.writer(csv_file)

for file in final_files:
    writer.writerow(file)    

