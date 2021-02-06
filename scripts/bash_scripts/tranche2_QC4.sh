vcftools --gzvcf snp.recalibrated.vcf.gz --remove-filtered-all --recode --out All_filters_pass --recode-INFO-all



cat All_filters_pass.recode.vcf| grep -v '#'| awk -F '\t' '{print $1, $2}' > pass_all_filters


vcftools --gzvcf /mnt/Genoma/drobles/ccastaneda/resources/gnomad.genomes.r3.0.sites.vcf.bgz --positions pass_all_filters --recode --recode-INFO-all --out gnomad_myregions_pass


vcftools --vcf gnomad_myregions_pass.recode.vcf --remove-filtered-all --recode --recode-INFO-all --out mycorrectSNPs


cat mycorrectSNPs.recode.vcf| grep -v '#'| awk -F '\t' '{print $1, $2}' > pass_gnomad_filters

vcftools --vcf ./All_filters_pass.recode.vcf --positions pass_gnomad_filters --recode  --recode-INFO-all --out filtered

plink --vcf ahora_si_ya_filtrado.recode.vcf --assoc fisher --pheno /mnt/Genoma/drobles/ccastaneda/gwas-reseq/archivos_fuente/pheno --out filtered --allow-no-sex --keep-allele-order --set-missing-var-ids @:#

cat filtered.assoc.fisher | tr -s ' ' ',' > filtered.fisher.csv

rm filtered.log
mm filtered.nosex 
mv filtered.assoc.fisher 
mv filtered.assoc.fisher assocs
mv filtered.fisher.csv assocs
rm All_filters_pass.recode.vcf
rm All_filters_pass.log 
rm gnomad*
rm mycorrectSNPs*
rm pass_gnomad_filters
rm cohort*
rm mdust*
rm pass*
rm super*
rm All*i
rm indel*
rm missing*
rm out*
rm snp*
mv plot* /mnt/Genoma/drobles/ccastaneda/gwas-reseq/gatk_plots

