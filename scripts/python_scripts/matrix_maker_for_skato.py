


import pandas as pd
import io
import os

vcf = 'filtered.recode.vcf'

with open(vcf, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
    vcf = pd.read_table(io.StringIO(str.join(os.linesep, lines)), dtype={'#CHROM':str, 'POS':str}, low_memory=False)
    vcf = vcf.rename(columns={'#CHROM': 'CHROM'})
    vcf = vcf.set_index(['CHROM', 'POS'])
    samples = vcf.iloc[:,7:].applymap(lambda s: str(s.split(':')[0]))


samples.replace(to_replace='1/1', value= '2', inplace=True, limit=None, regex=False, method='pad')
samples.replace(to_replace='0/0', value= '0', inplace=True, limit=None, regex=False, method='pad')
samples.replace(to_replace='0/1', value= '1', inplace=True, limit=None, regex=False, method='pad')
samples.replace(to_replace='1/0', value= '1', inplace=True, limit=None, regex=False, method='pad')
samples.replace(to_replace='1|0', value= '1', inplace=True, limit=None, regex=False, method='pad')
samples.replace(to_replace='0|0', value= '0', inplace=True, limit=None, regex=False, method='pad')
samples.replace(to_replace='0|1', value= '1', inplace=True, limit=None, regex=False, method='pad')
samples.replace(to_replace='1|1', value= '2', inplace=True, limit=None, regex=False, method='pad')


samples.to_csv('alleledosage.csv')
