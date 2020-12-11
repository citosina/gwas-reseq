
cd /mnt/Genoma/drobles/ccastaneda/gwas-reseq/lcr_filtered/regions

for file in $(find -name "*.recode.vcf");do plink --vcf $file --assoc fisher --pheno /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/pheno --out ${file/%.recode.vcf/} --allow-no-sex --keep-allele-order --set-missing-var-ids @:#;done


for file in $(find -name "*.assoc.fisher");do cat $file | tr -s ' ' ',' > ${file/.assoc.fisher/.assoc.csv};done



for file in $(find -name "*.nosex");do rm $file; done


