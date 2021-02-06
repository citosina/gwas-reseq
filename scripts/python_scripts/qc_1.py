
import pandas as pd
import numpy as np
import sys
import io 
import os

b1 = pd.read_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/samples_round1.txt', header = None)
b2 = pd.read_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/samples_round2.txt', header = None)
case = pd.read_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/case_samples.csv', header = None)
ctrl = pd.read_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/control_samples.csv', header = None)
ctrl_b1 = b1.merge(ctrl, on = 0)
ctrl_b2 = b2.merge(ctrl, on = 0)
case_b1 = b1.merge(case, on = 0)
case_b2 = b2.merge(case, on = 0)
ctrl_b1.to_csv('ctrl_b1', sep = "\t", index = False, header = None)
ctrl_b2.to_csv('ctrl_b2', sep = "\t", index = False, header = None)
case_b1.to_csv('case_b1', sep = "\t", index = False, header = None)
case_b2.to_csv('case_b2', sep = "\t", index = False, header = None)
print(ctrl_b1.info(), ctrl_b2.info(), case_b1.info(), case_b2.info())



