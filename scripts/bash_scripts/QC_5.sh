



source activate vep

/mnt/Genoma/drobles/ccastaneda/home/ensembl-vep/vep -i out.recode.vcf  --everything --cache -o filtered.vep  --output_format vcf --no_stats --force_overwrite

python vep_to_csv.py



