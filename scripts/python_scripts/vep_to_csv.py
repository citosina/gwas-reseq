



import pandas as pd
import io
import os


VEP = 'pasan.vep'



with open(VEP, 'r') as f:
    lines = [l for l in f if not l.startswith('##')]
    vep = pd.read_table(io.StringIO(str.join(os.linesep, lines)), dtype={'#CHROM':str, 'POS':str}, low_memory=False)
    vep = vep.rename(columns={'#CHROM': 'CHROM'})
    vep = vep.set_index(['CHROM', 'POS'])




vep = vep.reset_index()
vep['ID_x'] = vep.CHROM + ':' + vep.POS +':'+vep.REF+':'+vep.ALT
vep = vep.set_index(['CHROM', 'POS'])
IDx = vep.ID_x.to_frame()
variant_info = vep.iloc[:,:6]
variant_info = variant_info.drop(axis = 1, labels = ['ID', 'QUAL', 'FILTER'])
variant_info['samtools'] = variant_info.INFO.str.split('|').str.get(0)
variant_info['VEP'] = variant_info.INFO.str.split(';').str.get(-1)
variant_info = variant_info.reset_index().set_index(['CHROM', 'POS'])
VEP = variant_info.VEP.str.split(',', expand = True)
VEP = VEP.dropna(axis = 1, how = 'all')

stack = VEP.stack()


stack = stack.to_frame()


stack['Allele']= stack[0].str.split('|').str.get(0)
stack['Consequence']= stack[0].str.split('|').str.get(1)
stack['IMPACT']= stack[0].str.split('|').str.get(2)
stack['SYMBOL']= stack[0].str.split('|').str.get(3)
stack['Gene']= stack[0].str.split('|').str.get(4)
stack['Feature_type']= stack[0].str.split('|').str.get(5)
stack['Feature']= stack[0].str.split('|').str.get(6)
stack['BIOTYPE']= stack[0].str.split('|').str.get(7)
stack['EXON']= stack[0].str.split('|').str.get(8)
stack['INTRON']= stack[0].str.split('|').str.get(9)
stack['HGVSc']= stack[0].str.split('|').str.get(10)
stack['HGVSp']= stack[0].str.split('|').str.get(11)
stack['cDNA_position']= stack[0].str.split('|').str.get(12)
stack['CDS_position']= stack[0].str.split('|').str.get(13)
stack['Protein_position']= stack[0].str.split('|').str.get(14)
stack['Amino_acids']= stack[0].str.split('|').str.get(15)
stack['Codons']= stack[0].str.split('|').str.get(16)
stack['Existing_variation']= stack[0].str.split('|').str.get(17)
stack['DISTANCE']= stack[0].str.split('|').str.get(18)
stack['STRAND']= stack[0].str.split('|').str.get(19)
stack['FLAGS']= stack[0].str.split('|').str.get(20)
stack['VARIANT_CLASS']= stack[0].str.split('|').str.get(21)
stack['SYMBOL_SOURCE']= stack[0].str.split('|').str.get(22)
stack['HGNC_ID']= stack[0].str.split('|').str.get(23)
stack['CANONICAL']= stack[0].str.split('|').str.get(24)
stack['MANE']= stack[0].str.split('|').str.get(25)
stack['TSL']= stack[0].str.split('|').str.get(26)
stack['APPRIS']= stack[0].str.split('|').str.get(27)
stack['CCDS']= stack[0].str.split('|').str.get(28)
stack['ENSP']= stack[0].str.split('|').str.get(29)
stack['SWISSPROT']= stack[0].str.split('|').str.get(30)
stack['TREMBL']= stack[0].str.split('|').str.get(31)
stack['UNIPARC']= stack[0].str.split('|').str.get(32)
stack['GENE_PHENO']= stack[0].str.split('|').str.get(33)
stack['SIFT']= stack[0].str.split('|').str.get(34)
stack['PolyPhen']= stack[0].str.split('|').str.get(35)
stack['DOMAINS']= stack[0].str.split('|').str.get(36)
stack['miRNA']= stack[0].str.split('|').str.get(37)
stack['HGVS_OFFSET']= stack[0].str.split('|').str.get(38)
stack['AF']= stack[0].str.split('|').str.get(39)
stack['AFR_AF']= stack[0].str.split('|').str.get(40)
stack['AMR_AF']= stack[0].str.split('|').str.get(41)
stack['EAS_AF']= stack[0].str.split('|').str.get(42)
stack['EUR_AF']= stack[0].str.split('|').str.get(43)
stack['SAS_AF']= stack[0].str.split('|').str.get(44)
stack['AA_AF']= stack[0].str.split('|').str.get(45)
stack['EA_AF']= stack[0].str.split('|').str.get(46)
stack['gnomAD_AF']= stack[0].str.split('|').str.get(47)
stack['gnomAD_AFR_AF']= stack[0].str.split('|').str.get(48)
stack['gnomAD_AMR_AF']= stack[0].str.split('|').str.get(49)
stack['gnomAD_ASJ_AF']= stack[0].str.split('|').str.get(50)
stack['gnomAD_EAS_AF']= stack[0].str.split('|').str.get(51)
stack['gnomAD_FIN_AF']= stack[0].str.split('|').str.get(52)
stack['gnomAD_NFE_AF']= stack[0].str.split('|').str.get(53)
stack['gnomAD_OTH_AF']= stack[0].str.split('|').str.get(54)
stack['gnomAD_SAS_AF']= stack[0].str.split('|').str.get(55)
stack['MAX_AF']= stack[0].str.split('|').str.get(56)
stack['MAX_AF_POPS']= stack[0].str.split('|').str.get(57)
stack['CLIN_SIG']= stack[0].str.split('|').str.get(58)
stack['SOMATIC']= stack[0].str.split('|').str.get(59)
stack['PHENO']= stack[0].str.split('|').str.get(60)
stack['PUBMED']= stack[0].str.split('|').str.get(61)
stack['MOTIF_NAME']= stack[0].str.split('|').str.get(62)
stack['MOTIF_POS']= stack[0].str.split('|').str.get(63)
stack['HIGH_INF_POS']= stack[0].str.split('|').str.get(64)
stack['MOTIF_SCORE_CHANGE']= stack[0].str.split('|').str.get(65)


stack = stack.drop(axis = 1, labels = 0)

stack = stack.reset_index()

stack['SNP'] = stack.CHROM + ':' +stack.POS.astype(str)

stack = stack.set_index('CHROM', 'POS')
stack.to_csv('passvep_stack_vep.csv')

stack = stack.reset_index()
stack['BP'] = stack['POS']

chr1 = stack[(stack['CHROM'] == 'chr1')]
chr2 = stack[(stack['CHROM'] == 'chr2')]
chr3 = stack[(stack['CHROM'] == 'chr3')]
chr4 = stack[(stack['CHROM'] == 'chr4')]
chr5 = stack[(stack['CHROM'] == 'chr5')]
chr6 = stack[(stack['CHROM'] == 'chr6')]
chr7 = stack[(stack['CHROM'] == 'chr7')]
chr9 = stack[(stack['CHROM'] == 'chr9')]
chr10 = stack[(stack['CHROM'] == 'chr10')]
chr11 = stack[(stack['CHROM'] == 'chr11')]
chr12 = stack[(stack['CHROM'] == 'chr12')]
chr14 = stack[(stack['CHROM'] == 'chr14')]
chr15 = stack[(stack['CHROM'] == 'chr15')]
chr16 = stack[(stack['CHROM'] == 'chr16')]
chr17 = stack[(stack['CHROM'] == 'chr17')]
chr18 = stack[(stack['CHROM'] == 'chr18')]
chr20 = stack[(stack['CHROM'] == 'chr20')]
chr21 = stack[(stack['CHROM'] == 'chr21')]
chr22 = stack[(stack['CHROM'] == 'chr22')]
ARNT = chr1[(chr1['BP'].astype(int) > 150754900) & (chr1['BP'].astype(int) < 151073108)].rename(columns={'POS':'chr1'})
PARP1 = chr1[(chr1['BP'].astype(int) > 226258651) & (chr1['BP'].astype(int) < 226499883)].rename(columns={'POS':'chr1'})
CASP8 = chr2[(chr2['BP'].astype(int) > 201047085) & (chr2['BP'].astype(int) < 201724509)].rename(columns={'POS':'chr2'})
PAX3 = chr2[(chr2['BP'].astype(int) > 222095089) & (chr2['BP'].astype(int) < 222374945)].rename(columns={'POS':'chr2'})
BAP1 = chr3[(chr3['BP'].astype(int) > 52289996) & (chr3['BP'].astype(int) <  52520048)].rename(columns={'POS':'chr3'})
MITF = chr3[(chr3['BP'].astype(int) >  69636417) & (chr3['BP'].astype(int) <  70078609)].rename(columns={'POS':'chr3'})
SOD3 = chr4[(chr4['BP'].astype(int) > 24692289) & (chr4['BP'].astype(int) <  25055958)].rename(columns={'POS':'chr4'})
KIT = chr4[(chr4['BP'].astype(int) > 54548829) & (chr4['BP'].astype(int) <  54848706)].rename(columns={'POS':'chr4'})
GC = chr4[(chr4['BP'].astype(int) > 71703635) & (chr4['BP'].astype(int) <  71867003)].rename(columns={'POS':'chr4'})
NFKB1 = chr4[(chr4['BP'].astype(int) > 102394419) & (chr4['BP'].astype(int) < 102697771)].rename(columns={'POS':'chr4'})
TERT = chr5[(chr5['BP'].astype(int) > 1142197) & (chr5['BP'].astype(int) < 1447788)].rename(columns={'POS':'chr5'})
SLC45A2 = chr5[(chr5['BP'].astype(int) > 33840139) & (chr5['BP'].astype(int) < 34064027)].rename(columns={'POS':'chr5'})
IRF4 = chr6[(chr6['BP'].astype(int) > 281588) & (chr6['BP'].astype(int) < 512178)].rename(columns={'POS':'chr6'})
CYB5R4 = chr6[(chr6['BP'].astype(int) > 3948233) & (chr6['BP'].astype(int) <  84292890)].rename(columns={'POS':'chr6'})
IL6 = chr7[(chr7['BP'].astype(int) > 22617257) & (chr7['BP'].astype(int) <  22833613)].rename(columns={'POS':'chr7'})
MTAP = chr9[(chr9['BP'].astype(int) > 21665209) & (chr9['BP'].astype(int) <  22151817)].rename(columns={'POS':'chr9'})
NFKB2 = chr10[(chr10['BP'].astype(int) > 102289907) & (chr10['BP'].astype(int) < 102513171)].rename(columns={'POS':'chr10'})
CCND1 = chr11[(chr11['BP'].astype(int) > 69210414) & (chr11['BP'].astype(int) < 69803433)].rename(columns={'POS':'chr11'})
TYR = chr11[(chr11['BP'].astype(int) > 89097917) & (chr11['BP'].astype(int) <  89675699)].rename(columns={'POS':'chr11'})
ATM = chr11[(chr11['BP'].astype(int) > 108136364) & (chr11['BP'].astype(int) <  108468324)].rename(columns={'POS':'chr11'})
CYP27B1 = chr12[(chr12['BP'].astype(int) > 57666465 ) & (chr12['BP'].astype(int) < 57875869)].rename(columns={'POS':'chr12'})
TEP1 = chr14[(chr14['BP'].astype(int) > 20399357) & (chr14['BP'].astype(int) <  20550325)].rename(columns={'POS':'chr14'})
OCA2 = chr15[(chr15['BP'].astype(int) > 27670078) & (chr15['BP'].astype(int) < 28337214)].rename(columns={'POS':'chr15'})
GRIN2A = chr16[(chr16['BP'].astype(int) > 9696254) & (chr16['BP'].astype(int) <  10280553)].rename(columns={'POS':'chr16'})
SMG1 = chr16[(chr16['BP'].astype(int) >18749854 ) & (chr16['BP'].astype(int) < 19024774)].rename(columns={'POS':'chr16'})
FTO = chr16[(chr16['BP'].astype(int) > 53363904 ) & (chr16['BP'].astype(int) <  54384452)].rename(columns={'POS':'chr16'})
CDH1 = chr16[(chr16['BP'].astype(int) > 68528729) & (chr16['BP'].astype(int) < 68923897)].rename(columns={'POS':'chr16'})
MC1R = chr16[(chr16['BP'].astype(int) > 89806464) & (chr16['BP'].astype(int) < 90031424)].rename(columns={'POS':'chr16'})
CARD14 = chr17[(chr17['BP'].astype(int) > 80154433) & (chr17['BP'].astype(int) < 80508392)].rename(columns={'POS':'chr17'})
BCL2 = chr18[(chr18['BP'].astype(int) > 63015199) & (chr18['BP'].astype(int) < 63421327)].rename(columns={'POS':'chr18'})
ASIP = chr20[(chr20['BP'].astype(int) > 33709976) & (chr20['BP'].astype(int) < 34374270)].rename(columns={'POS':'chr20'})
TP53INP2 = chr20[(chr20['BP'].astype(int) > 34619476) & (chr20['BP'].astype(int) < 34933005)].rename(columns={'POS':'chr20'})
MX2 = chr21[(chr21['BP'].astype(int) > 41194545) & (chr21['BP'].astype(int) <  41512249)].rename(columns={'POS':'chr21'})
PLA2G6 = chr22[(chr22['BP'].astype(int) > 37870877) & (chr22['BP'].astype(int) <  38328731)].rename(columns={'POS':'chr22'})

ARNT.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/ARNT/ARNT_vep.csv')
PARP1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/PARP1/PARP1_vep.csv') 
CASP8.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/CASP8/CASP8_vep.csv') 
PAX3.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/PAX3/PAX3_vep.csv') 
BAP1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/BAP1/BAP1_vep.csv') 
MITF.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/MITF/MITF_vep.csv') 
SOD3.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/SOD3/SOD3_vep.csv') 
KIT.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/KIT/KIT_vep.csv') 
GC.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/GC/GC_vep.csv') 
NFKB1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/NFKB1/NFKB1_vep.csv') 
TERT.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/TERT/TERT_vep.csv') 
SLC45A2.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/SLC45A2/SLC45A2_vep.csv') 
IRF4.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/IRF4/IRF4_vep.csv') 
CYB5R4.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/CYB5R4/CYB5R4_vep.csv') 
IL6.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/IL6/IL6_vep.csv') 
MTAP.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/MTAP/MTAP_vep.csv') 
NFKB2.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/NFKB2/NFKB2_vep.csv') 
CCND1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/CCND1/CCND1_vep.csv') 
TYR.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/TYR/TYR_vep.csv') 
ATM.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/ATM/ATM_vep.csv') 
CYP27B1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/CYP27B1/CYP27B1_vep.csv')
TEP1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/TEP1/TEP1_vep.csv') 
OCA2.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/OCA2/OCA2_vep.csv') 
GRIN2A.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/GRIN2A/GRIN2A_vep.csv') 
SMG1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/SMG1/SMG1_vep.csv') 
FTO.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/FTO/FTO_vep.csv') 
CDH1.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/CDH1/CDH1_vep.csv') 
MC1R.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/MC1R/MC1R_vep.csv') 
CARD14.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/CARD14/CARD14_vep.csv') 
BCL2.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/BCL2/BCL2_vep.csv') 
ASIP.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/ASIP/ASIP_vep.csv') 
TP53INP2.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/TP53INP2/TP53INP2_vep.csv') 
MX2.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/MX2/MX2_vep.csv') 
PLA2G6.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/PLA2G6/PLA2G6_vep.csv') 
TERT.to_csv('/mnt/Genoma/drobles/ccastaneda/gwas-reseq/regions/TERT/TERT_vep.csv')


