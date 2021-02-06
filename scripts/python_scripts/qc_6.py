import pandas as pd
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import os
import math 
import numpy as np
import subprocess

vep = pd.read_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/passvep_stack_vep.csv', low_memory = False)
info = pd.read_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/filtered.fisher.csv')
gnomad = pd.read_csv('/mnt/Genoma/drobles/ccastaneda/resources/leeds_resources/gnomad_isin_reseq_af', delimiter = ',')


vep['SNP_reseq'] = vep.CHROM.astype(str) + ':' + vep.POS.astype(str)
vep = vep.set_index('SNP_reseq')
vep = vep.drop(axis = 1, labels = 'level_2')
print(vep.columns.values)
print(vep.head())


info['-logP'] = -(np.log10(info.P))
info['SNP_reseq'] = info.SNP
info = info.drop(axis = 1, labels = 'SNP')
info = info.set_index('SNP_reseq')
info = info.drop(axis = 1, labels = {'Unnamed: 0','Unnamed: 10'})
print(info.sort_values('P').head())


gnomad['SNP_reseq'] = gnomad['CHROM'] +':'+gnomad['POS_x'].astype(str)
gnomad['gnomad_AF'] = gnomad.AF_nfe
gnomad = gnomad.set_index('SNP_reseq')

everythin = vep.merge(info, right_index = True, left_index = True )
everything = everythin.merge(gnomad, right_index = True, left_index = True )
h = everything

sns.scatterplot(x="F_U", y='gnomad_AF', data=h, legend = False)
plt.savefig('F_Uvsgnomad_AF.tif')

h.to_csv('h.csv')

