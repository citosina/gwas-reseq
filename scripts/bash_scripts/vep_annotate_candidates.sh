



source activate vep
module load vcftools 

mkdir candidates

vcftools --vcf *.recode.vcf --positions candidates.bed --recode --recode-INFO-all --out ./candidates/candidates 

cd candidates

/mnt/Genoma/drobles/ccastaneda/home/ensembl-vep/vep -i candidates.recode.vcf --everything --plugin CADD --cache -o CADD --output_format tab --no_stats

/mnt/Genoma/drobles/ccastaneda/home/ensembl-vep/vep -i candidates.recode.vcf --everything --plugin Blosum62 --cache -o Blosum62 --output_format tab --no_stats

/mnt/Genoma/drobles/ccastaneda/home/ensembl-vep/vep -i candidates.recode.vcf --everything  --cache -o everything --output_format tab --no_stats

/mnt/Genoma/drobles/ccastaneda/home/ensembl-vep/vep -i candidates.recode.vcf --everything --plugin SameCodon --cache -o SameCodon --output_format tab --no_stats

/mnt/Genoma/drobles/ccastaneda/home/ensembl-vep/vep -i candidates.recode.vcf --plugin satMutMPRA,file=/home/ccastaneda/.vep/Plugins/satMutMPRA_GRCh38_ALL.gz,cols=ALL --cache -o satmutmpra --output_format tab --no_stats

cd ../

cat r1 r2 > s1
