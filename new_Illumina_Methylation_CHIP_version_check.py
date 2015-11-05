#http://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html

import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from numpy.random import randn
from scipy import stats
from pandas.tools.plotting import scatter_matrix
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import re
import os
import scipy.cluster.hierarchy as hier
import pylab as pl
import matplotlib.ticker as mticker
from Bio.pairwise2 import format_alignment
get_ipython().magic('matplotlib inline')


data_old = pd.read_csv('HumanMethylation450.csv', sep=',',  error_bad_lines=False, low_memory=False, skiprows=6)
data_old.columns=data_old.iloc[0] #Name the columns with the row[0]
data_old.dropna(subset=['CHR'], inplace=True)


data_old['Comb1'] = data_old['CHR'].map(str) + data_old['MAPINFO']


data_new = pd.read_csv('EPIC_manifest_file_1.csv', sep = ',', error_bad_lines = False, low_memory = False, skiprows=7)


data_new.dropna(subset=['CHR'],inplace=True)


def connect_two_columns(value):
    
    x, y = value
    new_name = str(x)+str(int(y))
    return(str(new_name))

data_new['Comb1'] = data_new[['CHR','MAPINFO']].apply(connect_two_columns, axis=1)



def remove_duplicates(values):
    output = []
    seen = set()
    for value in values:
        # If value has not been encountered yet,
        # ... add it to both list and set.
        if value not in seen:
            output.append(value)
            seen.add(value)
    return output


old_data_Comb1_List = data_old['Comb1'].tolist()

uniqueOld = remove_duplicates(old_data_Comb1_List)

new_data_Comb1_List = data_new['Comb1'].tolist()

uniqueNew = remove_duplicates(new_data_Comb1_List)

print("There are duplicates at OLD CHIP:", len(old_data_Comb1_List) - len(uniqueOld))
print("There are duplicates at NEW CHIP:", len(new_data_Comb1_List) - len(uniqueNew))


commonset_full = list(set(old_data_Comb1_List) & set(new_data_Comb1_List))

dif_set_old = list(set(old_data_Comb1_List) - set(commonset_full))

dif_set_new = list(set(new_data_Comb1_List) - set(commonset_full))

CommonDF = data_new[data_new['Comb1'].isin(commonset_full)]
CommonDF = CommonDF.ix[1:]

difOldDF = data_old[data_old['Comb1'].isin(dif_set_old)]
difOldDF = difOldDF.ix[1:]

difNewDF = data_new[data_new['Comb1'].isin(dif_set_new)]
difNewDF = difNewDF.ix[1:]

difOldDF.groupby('Regulatory_Feature_Group').size()

difOldDF.groupby('Enhancer').size()

difOldDF.groupby('Random_Loci').size()

difOldDF.groupby('Methyl27_Loci').size()

difOldDF.groupby('DMR').size()

difOldDF_genes = difOldDF.dropna(subset = ['UCSC_RefGene_Name'])
print("Number of entries with Gene reference:", len(difOldDF_genes))

commonN = len(commonset_full)
uniqueOld = len(dif_set_old)
uniqueNew = len(dif_set_new)

print("Old CHIP size:", len(old_data_Comb1_List))
print("New CHIP size:", len(new_data_Comb1_List))
print("Common Set:", commonN)
print("Unique at OLD ChIP:", uniqueOld)
print("Unique at New CHIP:", uniqueNew)

def save_DF_to_csv(df, fileName, folderName ):
    
    make_directory(folderName)
    path = folderName + '/' + fileName + '.csv'
    df.to_csv(os.path.join(path), sep = ',')

#---------------------------------------------------
def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
#---------------------------------------------------
folderName = 'Results'

fileName1 = 'Common_Set'
df1 = CommonDF

fileName2 = 'UniqueOldChip'
df2 = difOldDF

fileName3 = 'UniqueNewChip'
df3 = difNewDF

save_DF_to_csv(df1, fileName1 , folderName)
save_DF_to_csv(df2, fileName2 , folderName)
save_DF_to_csv(df3, fileName3 , folderName)
