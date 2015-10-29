
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
small_data_old = data_old[['IlmnID','CHR', 'MAPINFO','Forward_Sequence', 'Strand', 'UCSC_RefGene_Name']]

#combine columns
small_data_old['Comb1'] = small_data_old['CHR'].map(str) + small_data_old['MAPINFO']

small_data_old['CombFull'] = small_data_old['Comb1'].map(str) + small_data_old['Strand']


small_data_old_genes = small_data_old.dropna(subset = ['UCSC_RefGene_Name'])
small_data_old_others = small_data_old[small_data_old['UCSC_RefGene_Name'].isnull()]



data_new = pd.read_csv('HumanMethylation550.csv', sep = ',', error_bad_lines = False, low_memory = False, skiprows=6)
data_new.columns=data_new.iloc[0]
data_new.dropna(subset=['CHR'],inplace=True)
small_data_new = data_new[['IlmnID','CHR', 'MAPINFO','Forward_Sequence', 'Strand', 'UCSC_RefGene_Name']]



small_data_new['Comb1'] = small_data_new['CHR'].map(str) + small_data_new['MAPINFO']
small_data_new['CombFull'] = small_data_new['Comb1'].map(str) + small_data_new['Strand']
small_data_new_genes = small_data_new.dropna(subset = ['UCSC_RefGene_Name'])
small_data_new_others = small_data_new[small_data_old['UCSC_RefGene_Name'].isnull()]


small_data_new_others.shape

old_genes_Combfull = small_data_old_genes['CombFull'].tolist()

new_genes_Combfull = small_data_new_genes['CombFull'].tolist()

Common_gene_CombFull = list(set(old_genes_Combfull) & set(new_genes_Combfull))

Common_new_genes_Combfull_DF_1 = small_data_new_genes[small_data_new_genes['CombFull'].isin(Common_gene_CombFull)]

Common_old_genes_Combfull_DF_1 = small_data_old_genes[small_data_old_genes['CombFull'].isin(Common_gene_CombFull)]

old_gene_Name_select2 = small_data_old_genes['UCSC_RefGene_Name'].tolist()


new_gene_Name_select2 = small_data_new_genes['UCSC_RefGene_Name'].tolist()

Common_select2 = list(set(old_genes_Combfull) & set(new_genes_Combfull))


if len(Common_select2) == len(Common_gene_CombFull):
    print("No difference between mapping and gene names %5d - %5d" %(len(Common_select2), len(Common_gene_CombFull)))
else:
    print ("There is a difference between mapping and gene names. \n Difference = %5d" % (len(Common_gene_CombFull) - len(Common_select2)))


dif_old_gene_CombFull = list(set(old_genes_Combfull) - set(new_genes_Combfull))


dif_new_gene_CombFull = list(set(new_genes_Combfull) - set(old_genes_Combfull))


dif_old_gene = small_data_old_genes[small_data_old_genes['CombFull'].isin(dif_old_gene_CombFull)]


dif_new_gene = small_data_new_genes[small_data_new_genes['CombFull'].isin(dif_new_gene_CombFull)]


chrList = list(range(1,23))


dif_old_mapinfo = dif_old_gene['MAPINFO'].tolist()

dif_new_mapinfo = dif_new_gene['MAPINFO'].tolist()


def plotting_line(x1,y1,x2,y2, ChrName):
    
    fullList = list(x1)
    fullList.extend(list(x2))
    fullList.sort()
    fmax = max(fullList)
    fmin = min(fullList)
    
    fig = plt.figure()
    ax1 = plt.subplot2grid((1,1),(0,0))
    
    plt.scatter(x1,y1, facecolor='g', label = 'New')
    plt.scatter(x2,y2, label = 'Old')
    
    plt.xlabel('Chr MAPINFO')
    plt.ylabel('1: new, 2:old')
    
    title = 'Cromosome: ' + str(ChrName)
    plt.title(title)
    
    ax1.yaxis.label.set_color('m')
    ax1.xaxis.label.set_color('c')
    ax1.set_yticks([1,2])
    ax1.xaxis.set_major_locator(mticker.MaxNLocator(30))
    #ax1.set_xticks(fullList)
    
    for label in ax1.xaxis.get_ticklabels():
        label.set_rotation(45)
    
    
    ax1.grid(True)
    plt.subplots_adjust(left=.09, bottom = .16, right = .94, top =.94, wspace = .20, hspace = .2)
    plt.legend()
    plt.show()


list_old = list([2]*len(dif_old_mapinfo))
list_new = list([1]*len(dif_new_mapinfo))

chromList = list(range(1,23))

extra = ['Y','X']

chromList.extend(extra)

for chrom in chromList:
    
    nx1 = 0
    nx2 = 0
    tempDF1 = []
    tempDF2 = []
    
    tempDF1 = dif_new_gene[dif_new_gene['CHR'] == str(chrom)]
    y1 = tempDF1['MAPINFO'].tolist()
    
    nx1 = len(y1)
    x1 = list([1]*len(y1))
    
    tempDF2 = dif_old_gene[dif_old_gene['CHR'] == str(chrom)]
    y2 = tempDF2['MAPINFO'].tolist()
    
    nx2 = len(y2)
    x2 = list([2]*len(y2))
    
    print(chrom, nx1, nx2)
    
    if nx1 > 0 or nx2 > 0:
        plotting_line(y1, x1, y2, x2, chrom)
        

plotting_line(dif_old_mapinfo, list_old, dif_new_mapinfo, list_new, 1)
