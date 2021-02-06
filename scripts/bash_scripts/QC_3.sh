


module load bedtools
module load plink2/1.90beta3

vcftools --vcf cohort.recode.vcf --positions pass_almost_all_filters --out pass_almost_all_filters --recode --recode-INFO-all     
bedtools coverage -a pass_almost_all_filters.recode.vcf -b /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/mdust38.bed > mdust_overlap
cat mdust_overlap | awk -F '\t' '{print $1, $2, $4, $5,  $(NF-3), $(NF-2), $(NF-1), $NF}' > mdust_overlap_val

python /mnt/Genoma/drobles/ccastaneda/gwas-reseq/python_scripts/qc_4.py
echo "qc_4 ok \n output mdust_vars.csv"

awk -F ',' '{print $1, $2}' mdust_vars.csv > mdust_overlap_poss
vcftools --vcf pass_almost_all_filters.recode.vcf --exclude-positions mdust_overlap_poss --recode-INFO-all --recode --out super_filtered
cat super_filtered.recode.vcf| sed s/SC_//g > super_filtered_nwnms.vcf


plink --vcf super_filtered_nwnms.vcf --recode --out cohort --make-bed --pheno /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/pheno --allow-no-sex --no-fid --no-parents --set-missing-var-ids @:# --keep-allele-order 
plink --bfile cohort --missing --out cohort --pheno /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/pheno --allow-no-sex --set-missing-var-ids @:# --keep-allele-order
plink --bfile cohort --het --out cohort --pheno /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/pheno --allow-no-sex --set-missing-var-ids @:# --keep-allele-order
plink --bfile cohort --test-missing midp --out cohort --allow-no-sex --set-missing-var-ids @:#:\$1:\$2 --keep-allele-order 


python /mnt/Genoma/drobles/ccastaneda/gwas-reseq/python_scripts/qc_5.py
echo "qc_5 ok output missing_no_pass"

vcftools --vcf super_filtered_nwnms.vcf --exclude-positions missing_no_pass --recode --recode-INFO-all

