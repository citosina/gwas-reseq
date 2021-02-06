

module load vcftools
module load plink2
module load bedtools
module load gatk 

vcftools --gzvcf /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/all.vcf.gz --keep /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/cohort --out cohort --recode-INFO-all --recode
vcftools --vcf cohort.recode.vcf --keep /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/samples_round1.txt --recode-INFO-all --recode --out round1
vcftools --vcf cohort.recode.vcf --keep /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/samples_round2.txt --recode-INFO-all  --recode --out round2

python /mnt/Genoma/drobles/ccastaneda/gwas-reseq/python_scripts/qc_1.py
echo "qc_1.py oki \n outputs ctrl_b1 ctrl_b2 case_b1 case_b2"


