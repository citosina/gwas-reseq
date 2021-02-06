


vcftools --vcf ctrl_b1_hwe.recode.vcf --exclude-positions ctrl_b1_gq_avgPOS.tsv --recode --recode-INFO-all --out ctrl_b1_after_gqavg
vcftools --vcf ctrl_b2_hwe.recode.vcf --exclude-positions ctrl_b2_gq_avgPOS.tsv --recode --recode-INFO-all --out ctrl_b2_after_gqavg
vcftools --vcf case_b1.recode.vcf --exclude-positions case_b1_gq_avgPOS.tsv --recode --recode-INFO-all --out case_b1_after_gqavg
vcftools --vcf case_b2.recode.vcf --exclude-positions case_b2_gq_avgPOS.tsv --recode --recode-INFO-all --out  case_b2_after_gqavg


vcftools --vcf ctrl_b1_after_gqavg.recode.vcf --max-missing 0.05 --recode-INFO-all --recode --out ctrl_b1_POSfiltered
vcftools --vcf ctrl_b2_after_gqavg.recode.vcf --max-missing 0.05 --recode-INFO-all --recode --out ctrl_b2_POSfiltered
vcftools --vcf case_b1_after_gqavg.recode.vcf --max-missing 0.05 --recode-INFO-all --recode --out case_b1_POSfiltered
vcftools --vcf case_b2_after_gqavg.recode.vcf --max-missing 0.05 --recode-INFO-all --recode --out case_b2_POSfiltered


grep -v '^\#' ctrl_b1_POSfiltered.recode.vcf | awk -F'\t' '{ print $1, '\t', $2}' > ctrl_b1_POS.csv
grep -v '^\#' ctrl_b2_POSfiltered.recode.vcf | awk -F'\t' '{ print $1, '\t', $2}' > ctrl_b2_POS.csv
grep -v '^\#' case_b1_POSfiltered.recode.vcf | awk -F'\t' '{ print $1, '\t', $2}' > case_b1_POS.csv
grep -v '^\#' case_b2_POSfiltered.recode.vcf | awk -F'\t' '{ print $1, '\t', $2}' > case_b2_POS.csv
grep -v '^\#' cohort.recode.vcf | awk -F'\t' '{ print $1, '\t', $2}' > cohort_POS.csv


python /mnt/Genoma/drobles/ccastaneda/gwas-reseq/python_scripts/qc_3.py
echo "qc_3 ok \n output pass_almost_all_filters "

