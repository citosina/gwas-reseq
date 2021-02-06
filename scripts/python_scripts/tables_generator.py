


import os
import io
import pandas as pd
import numpy as np
import subprocess
import matplotlib.pyplot as plt

region = str(input("enter the region you want to figure:"))
wd = '/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/'+region
os.chdir(wd)

r1 = pd.read_csv('r1')
r2 = pd.read_csv('r2')

c1 = r1[(r1.Protein_position.astype(str) != 'nan' )]
r1 = r1.drop_duplicates('Existing_variation')

c2 = r2[(r2.Protein_position.astype(str) != 'nan' )]
r2 = r2.drop_duplicates('Existing_variation')

s = r1.merge(r2, how = 'outer')

c = c1.merge(c2, how = 'outer')

c['Position'] = c.GRCh38.str.split('-').str.get(0)

c = c.set_index('Position')

c = c.drop_duplicates('Existing_variation')

c['Protein_position'] =c.Protein_position.apply(lambda x: int(x))

fisher = pd.read_csv(region +'_fisher.csv').set_index('SNP')

col_list = ["SNP", "Protein_position", 'Amino_acids']
vep = pd.read_csv(region +'_vep.csv', usecols=col_list, low_memory = False).set_index('SNP')

vep = vep.dropna(axis = 0)

vep = vep.drop_duplicates()

g = pd.DataFrame()

g['Position'] = s.GRCh38.str.split('-').str.get(0)
g['Base change'] = s.ref_alt.apply(lambda x: str(x).replace('_', '/'))

g = g.fillna(np.NaN)

g['Base change'] = np.where(g['Base change'] == 'nan', s['ref_alt_x'].apply(lambda x: str(x).replace('_', '/')), g['Base change'])

g = g.set_index('Position')

g['AF in cases'] = fisher['F_A']

g['AF in controls'] = fisher['F_U']

g['P value'] = fisher['P']

g['OR'] = fisher['OR']

g = g.sort_values('P value')

g['Consequence'] = c.Amino_acids.str.split('/').str.get(0) + c.Protein_position.astype(str) + c.Amino_acids.str.split('/').str.get(1)

to_condition = pd.read_csv('to_condition', header = None, names = {'Position'}).set_index('Position')

to_condition['Base change'] = fisher.A1 +'/'+ fisher.A2

to_condition['AF in cases'] = fisher.F_A
to_condition['AF in controls'] = fisher.F_U

to_condition['P value'] = fisher.P
to_condition['OR'] = fisher.OR

vep['Consequence'] = vep.Amino_acids.str.split('/').str.get(0) + vep.Protein_position.astype(str) + vep.Amino_acids.str.split('/').str.get(1)

vep = vep.reset_index().rename(columns = {'SNP': 'Position'}).set_index('Position').drop\
(axis = 1, columns = {'Protein_position', 'Amino_acids'})

to_condition = to_condition.merge(vep, right_index = True, left_index = True)


s['Position'] = s.GRCh38.str.split('-').str.get(0)
s = s.set_index('Position')

g['P value 1st'] = s.P_cond1

g['P value 2nd'] = s.P_cond2

g['Region'] = s.Signal.str.replace('s', 'r')

add1 = pd.read_csv('./plink/1ADD.csv', header = None, names = [0, 'chr', 'Position', 'BP', 'A', 'test',\
                                                               'N', 'OR', 'STAT', 'P'])

add1 = add1.set_index('Position')

add2 = pd.read_csv('./plink/2ADD.csv', header = None, names = [0, 'chr', 'Position', 'BP', 'A', 'test',\
                                                               'N', 'OR', 'STAT', 'P'])

add2 = add2.set_index('Position')

to_condition['P value 1st'] = add1.P
to_condition['P value 2nd'] = add2.P

to_condition['Region'] = (range(1, 1 + len(to_condition)))
to_condition['Region'] = 'r' + to_condition.Region.astype(str) + ' lead'

to_condition

frames = [g, to_condition]
result = pd.concat(frames)

result = result.sort_values('P value').loc[~result.index.duplicated(keep='last')]

result.to_csv('table_for_paper.csv')


add1['-logP'] = -(np.log10(add1.P))
add2['-logP'] = -(np.log10(add2.P))

fig, ax = plt.subplots(3,1, sharey=True, sharex=True)

ax[0].scatter(fisher['BP'], fisher['-logP'])
ax[1].scatter(add1['BP'], add1['-logP'])
ax[2].scatter(add2['BP'], add2['-logP'])
plt.show()
plt.savefig('regresions.tiff')



