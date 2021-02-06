

import os
import io
import pandas as pd
import numpy as np

exons = pd.read_csv('../archivos_fuente/all_exons_design', header = None, delimiter = '\s+')

exons = exons.rename(columns = {0:'chr', 1:'start', 2:'end', 3:'gene'})

genes = pd.read_csv('../archivos_fuente/genes_start_end_strand', delimiter = '\t', header = None, names = ['Gene_id', 'gene','tss', 'gene_end', 'Strand'])

genes = genes.set_index('gene')

exons = exons.set_index('gene')

c = genes.merge(exons, right_index = True, left_index = True)

c['start_factor'] = c.Strand.apply(lambda x: -5*x)
c['end_factor'] = c.Strand.apply(lambda x: 5*x)

c.start = c.start + c.start_factor

c['promoter_factor'] =c.Strand.apply(lambda x: -200*x)
c['promoter'] = c['tss'] + c.promoter_factor

c.to_csv('skato-subsets1.csv')
