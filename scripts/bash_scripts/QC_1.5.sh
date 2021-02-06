

python /mnt/Genoma/drobles/ccastaneda/gwas-reseq/python_scripts/qc_1.py
echo "qc_1.py oki \n outputs ctrl_b1 ctrl_b2 case_b1 case_b2"


vcftools --vcf round1.recode.vcf --keep ctrl_b1 --recode --recode-INFO-all --out ctrl_b1
vcftools --vcf round2.recode.vcf --keep ctrl_b2 --recode --recode-INFO-all --out ctrl_b2
vcftools --vcf round1.recode.vcf --keep case_b1 --recode --recode-INFO-all --out case_b1
vcftools --vcf round2.recode.vcf --keep case_b2 --recode --recode-INFO-all --out case_b2



vcftools --vcf ctrl_b1.recode.vcf --minGQ 20 --recode --recode-INFO-all --out ctrl_b1_GQ
vcftools --vcf ctrl_b2.recode.vcf --minGQ 20 --recode --recode-INFO-all --out ctrl_b2_GQ
vcftools --vcf case_b1.recode.vcf --minGQ 20 --recode --recode-INFO-all --out ctrl_b1_GQ
vcftools --vcf case_b2.recode.vcf --minGQ 20 --recode --recode-INFO-all --out ctrl_b2_GQ




vcftools --vcf ctrl_b1.recode.vcf --hardy --recode-INFO-all --out ctrl_b1_hwe
vcftools --vcf ctrl_b2.recode.vcf --hardy --recode-INFO-all --out ctrl_b2_hwe
vcftools --vcf case_b1.recode.vcf --hardy --recode-INFO-all --out case_b1_hwe
vcftools --vcf case_b2.recode.vcf --hardy --recode-INFO-all --out case_b2_hwe



vcftools --vcf ctrl_b1.recode.vcf --hwe 0.001 --recode --recode-INFO-all --out ctrl_b1_hwe
vcftools --vcf ctrl_b2.recode.vcf --hwe 0.001 --recode --recode-INFO-all --out ctrl_b2_hwe



python /mnt/Genoma/drobles/ccastaneda/gwas-reseq/python_scripts/qc_2.py
echo "qc_2.py ok \n output ctrl_b1_gq_avgPOS.tsv ctrl_b2_gq_avgPOS.tsv  case_b1_gq_avgPOS.tsv case_b2_gq_avgPOS.tsv"

