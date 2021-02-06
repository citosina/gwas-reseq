


import pandas as pd
import numpy as np
a = pd.read_csv( 'pasan.fisher.csv', low_memory = False, index_col = None)
a['CHR'] = a.CHR.apply(lambda s: int(s))
a['ID_x'] = a.SNP 
a = a.drop(axis = 1, labels = 'Unnamed: 0')
a['POS'] = a.BP
a['CHROM'] = 'chr' + a.CHR.astype(str)
a = a.set_index('SNP')
a['-logP'] = -np.log10(a.P)



chr1 = a[(a['CHROM'] == 'chr1')]
chr2 = a[(a['CHROM'] == 'chr2')]
chr3 = a[(a['CHROM'] == 'chr3')]
chr4 = a[(a['CHROM'] == 'chr4')]
chr5 = a[(a['CHROM'] == 'chr5')]
chr6 = a[(a['CHROM'] == 'chr6')]
chr7 = a[(a['CHROM'] == 'chr7')]
chr9 = a[(a['CHROM'] == 'chr9')]
chr10 = a[(a['CHROM'] == 'chr10')]
chr11 = a[(a['CHROM'] == 'chr11')]
chr12 = a[(a['CHROM'] == 'chr12')]
chr14 = a[(a['CHROM'] == 'chr14')]
chr15 = a[(a['CHROM'] == 'chr15')]
chr16 = a[(a['CHROM'] == 'chr16')]
chr17 = a[(a['CHROM'] == 'chr17')]
chr18 = a[(a['CHROM'] == 'chr18')]
chr20 = a[(a['CHROM'] == 'chr20')]
chr21 = a[(a['CHROM'] == 'chr21')]
chr22 = a[(a['CHROM'] == 'chr22')]
ARNT = chr1[(chr1['BP'].astype(int) > 150754900) & (chr1['BP'].astype(int) < 151073108)].sort_values('P').rename(columns={'POS':'chr1'})
PARP1 = chr1[(chr1['BP'].astype(int) > 226258651) & (chr1['BP'].astype(int) < 226499883)].sort_values('P').rename(columns={'POS':'chr1'})
CASP8 = chr2[(chr2['BP'].astype(int) > 201047085) & (chr2['BP'].astype(int) < 201724509)].sort_values('P').rename(columns={'POS':'chr2'})
PAX3 = chr2[(chr2['BP'].astype(int) > 222095089) & (chr2['BP'].astype(int) < 222374945)].sort_values('P').rename(columns={'POS':'chr2'})
BAP1 = chr3[(chr3['BP'].astype(int) > 52289996) & (chr3['BP'].astype(int) <  52520048)].sort_values('P').rename(columns={'POS':'chr3'})
MITF = chr3[(chr3['BP'].astype(int) >  69636417) & (chr3['BP'].astype(int) <  70078609)].sort_values('P').rename(columns={'POS':'chr3'})
SOD3 = chr4[(chr4['BP'].astype(int) > 24692289) & (chr4['BP'].astype(int) <  25055958)].sort_values('P').rename(columns={'POS':'chr4'})
KIT = chr4[(chr4['BP'].astype(int) > 54548829) & (chr4['BP'].astype(int) <  54848706)].sort_values('P').rename(columns={'POS':'chr4'})
GC = chr4[(chr4['BP'].astype(int) > 71703635) & (chr4['BP'].astype(int) <  71867003)].sort_values('P').rename(columns={'POS':'chr4'})
NFKB1 = chr4[(chr4['BP'].astype(int) > 102394419) & (chr4['BP'].astype(int) < 102697771)].sort_values('P').rename(columns={'POS':'chr4'})
TERT = chr5[(chr5['BP'].astype(int) > 1142197) & (chr5['BP'].astype(int) < 1447788)].sort_values('P').rename(columns={'POS':'chr5'})
SLC45A2 = chr5[(chr5['BP'].astype(int) > 33840139) & (chr5['BP'].astype(int) < 34064027)].sort_values('P').rename(columns={'POS':'chr5'})
IRF4 = chr6[(chr6['BP'].astype(int) > 281588) & (chr6['BP'].astype(int) < 512178)].sort_values('P').rename(columns={'POS':'chr6'})
CYB5R4 = chr6[(chr6['BP'].astype(int) > 3948233) & (chr6['BP'].astype(int) <  84292890)].sort_values('P').rename(columns={'POS':'chr6'})
IL6 = chr7[(chr7['BP'].astype(int) > 22617257) & (chr7['BP'].astype(int) <  22833613)].sort_values('P').rename(columns={'POS':'chr7'})
MTAP = chr9[(chr9['BP'].astype(int) > 21665209) & (chr9['BP'].astype(int) <  22151817)].sort_values('P').rename(columns={'POS':'chr9'})
NFKB2 = chr10[(chr10['BP'].astype(int) > 102289907) & (chr10['BP'].astype(int) < 102513171)].sort_values('P').rename(columns={'POS':'chr10'})
CCND1 = chr11[(chr11['BP'].astype(int) > 69210414) & (chr11['BP'].astype(int) < 69803433)].sort_values('P').rename(columns={'POS':'chr11'})
TYR = chr11[(chr11['BP'].astype(int) > 89097917) & (chr11['BP'].astype(int) <  89675699)].sort_values('P').rename(columns={'POS':'chr11'})
ATM = chr11[(chr11['BP'].astype(int) > 108136364) & (chr11['BP'].astype(int) <  108468324)].sort_values('P').rename(columns={'POS':'chr11'})
CYP27B1 = chr12[(chr12['BP'].astype(int) > 57666465 ) & (chr12['BP'].astype(int) < 57875869)].sort_values('P').rename(columns={'POS':'chr12'})
TEP1 = chr14[(chr14['BP'].astype(int) > 20399357) & (chr14['BP'].astype(int) <  20550325)].sort_values('P').rename(columns={'POS':'chr14'})
OCA2 = chr15[(chr15['BP'].astype(int) > 27670078) & (chr15['BP'].astype(int) < 28337214)].sort_values('P').rename(columns={'POS':'chr15'})
GRIN2A = chr16[(chr16['BP'].astype(int) > 9696254) & (chr16['BP'].astype(int) <  10280553)].sort_values('P').rename(columns={'POS':'chr16'})
SMG1 = chr16[(chr16['BP'].astype(int) >18749854 ) & (chr16['BP'].astype(int) < 19024774)].sort_values('P').rename(columns={'POS':'chr16'})
FTO = chr16[(chr16['BP'].astype(int) > 53363904 ) & (chr16['BP'].astype(int) <  54384452)].sort_values('P').rename(columns={'POS':'chr16'})
CDH1 = chr16[(chr16['BP'].astype(int) > 68528729) & (chr16['BP'].astype(int) < 68923897)].sort_values('P').rename(columns={'POS':'chr16'})
MC1R = chr16[(chr16['BP'].astype(int) > 89806464) & (chr16['BP'].astype(int) < 90031424)].sort_values('P').rename(columns={'POS':'chr16'})
CARD14 = chr17[(chr17['BP'].astype(int) > 80154433) & (chr17['BP'].astype(int) < 80508392)].sort_values('P').rename(columns={'POS':'chr17'})
BCL2 = chr18[(chr18['BP'].astype(int) > 63015199) & (chr18['BP'].astype(int) < 63421327)].sort_values('P').rename(columns={'POS':'chr18'})
ASIP = chr20[(chr20['BP'].astype(int) > 33709976) & (chr20['BP'].astype(int) < 34374270)].sort_values('P').rename(columns={'POS':'chr20'})
TP53INP2 = chr20[(chr20['BP'].astype(int) > 34619476) & (chr20['BP'].astype(int) < 34933005)].sort_values('P').rename(columns={'POS':'chr20'})
MX2 = chr21[(chr21['BP'].astype(int) > 41194545) & (chr21['BP'].astype(int) <  41512249)].sort_values('P').rename(columns={'POS':'chr21'})
PLA2G6 = chr22[(chr22['BP'].astype(int) > 37870877) & (chr22['BP'].astype(int) <  38328731)].sort_values('P').rename(columns={'POS':'chr22'})


ARNT.sort_values('P').to_csv('./regions/ARNT/ARNT_fisher.csv')
PARP1.sort_values('P').to_csv('./regions/PARP1/PARP1_fisher.csv') 
CASP8.sort_values('P').to_csv('./regions/CASP8/CASP8_fisher.csv') 
PAX3.sort_values('P').to_csv('./regions/PAX3/PAX3_fisher.csv') 
BAP1.sort_values('P').to_csv('./regions/BAP1/BAP1_fisher.csv') 
MITF.sort_values('P').to_csv('./regions/MITF/MITF_fisher.csv') 
SOD3.sort_values('P').to_csv('./regions/SOD3/SOD3_fisher.csv') 
KIT.sort_values('P').to_csv('./regions/KIT/KIT_fisher.csv') 
GC.sort_values('P').to_csv('./regions/GC/GC_fisher.csv') 
NFKB1.sort_values('P').to_csv('./regions/NFKB1/NFKB1_fisher.csv') 
TERT.sort_values('P').to_csv('./regions/TERT/TERT_fisher.csv') 
SLC45A2.sort_values('P').to_csv('./regions/SLC45A2/SLC45A2_fisher.csv') 
IRF4.sort_values('P').to_csv('./regions/IRF4/IRF4_fisher.csv') 
CYB5R4.sort_values('P').to_csv('./regions/CYB5R4/CYB5R4_fisher.csv') 
IL6.sort_values('P').to_csv('./regions/IL6/IL6_fisher.csv') 
MTAP.sort_values('P').to_csv('./regions/MTAP/MTAP_fisher.csv') 
NFKB2.sort_values('P').to_csv('./regions/NFKB2/NFKB2_fisher.csv') 
CCND1.sort_values('P').to_csv('./regions/CCND1/CCND1_fisher.csv') 
TYR.sort_values('P').to_csv('./regions/TYR/TYR_fisher.csv') 
ATM.sort_values('P').to_csv('./regions/ATM/ATM_fisher.csv') 
CYP27B1.sort_values('P').to_csv('./regions/CYP27B1/CYP27B1_fisher.csv')
TEP1.sort_values('P').to_csv('./regions/TEP1/TEP1_fisher.csv') 
OCA2.sort_values('P').to_csv('./regions/OCA2/OCA2_fisher.csv') 
GRIN2A.sort_values('P').to_csv('./regions/GRIN2A/GRIN2A_fisher.csv') 
SMG1.sort_values('P').to_csv('./regions/SMG1/SMG1_fisher.csv') 
FTO.sort_values('P').to_csv('./regions/FTO/FTO_fisher.csv') 
CDH1.sort_values('P').to_csv('./regions/CDH1/CDH1_fisher.csv') 
MC1R.sort_values('P').to_csv('./regions/MC1R/MC1R_fisher.csv') 
CARD14.sort_values('P').to_csv('./regions/CARD14/CARD14_fisher.csv') 
BCL2.sort_values('P').to_csv('./regions/BCL2/BCL2_fisher.csv') 
ASIP.sort_values('P').to_csv('./regions/ASIP/ASIP_fisher.csv') 
TP53INP2.sort_values('P').to_csv('./regions/TP53INP2/TP53INP2_fisher.csv') 
MX2.sort_values('P').to_csv('./regions/MX2/MX2_fisher.csv') 
PLA2G6.sort_values('P').to_csv('./regions/PLA2G6/PLA2G6_fisher.csv') 



