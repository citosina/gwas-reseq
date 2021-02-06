









import os
import io
import pandas as pd
import numpy as np
import subprocess

region = str(input("enter the region you want to regress: "))
chromosome = str(input("enter the chromosome that region is at : "))

wd = '/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/'+region+'/candidates'

os.chdir(wd)


candidates = 'candidates.recode.vcf'
with open (candidates, 'r') as f:
    lines = [l for l in f if not l.startswith('##')] 
    variants = pd.read_table(io.StringIO(str.join(os.linesep, lines)), low_memory=False)
    variants = variants.rename(columns={'#CHROM': 'CHROM'})

variants['ALT'] = variants.ALT.str.replace('*', 'DEL')

variants['SNP'] = variants.CHROM + ':' + variants.POS.astype(str) +':' + variants.REF + ':' + variants.ALT
variants = variants.set_index('SNP')

who = variants.iloc[:,9:].applymap(lambda s: str(s.split(':')[0]))

who = who.reset_index()
who['CHROM'] = who.SNP.str.split(':').str.get(0)
who['POS'] = who.SNP.str.split(':').str.get(1).astype(int)
who = who.set_index(['CHROM', 'POS'])

a = who

a['ALT'] = a.SNP.str.split(':').str.get(-1)

a = a.replace(to_replace = np.NaN, value = 'no_rs')

a = a.replace(to_replace = '0/0', value=0)
a = a.replace(to_replace = '0/1', value=1)
a = a.replace(to_replace = '1/0', value=1)
a = a.replace(to_replace = '0|0', value=0)
a = a.replace(to_replace = '0|1', value=1)
a = a.replace(to_replace = '1|0', value=1)

a = a.reset_index().drop(axis = 1, labels = ['CHROM', 'POS']).set_index('SNP')

s1 = a.apply(lambda x: ','.join(x.index[x == 1]), axis=1)

s1 = s1.to_frame()

s1 = s1.rename(columns = {0:'carriers'})
print('These are the candidate carriers:')
print(s1)
s1.to_csv('candidates_carriers.tsv', sep = '\t')

