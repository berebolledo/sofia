#! /bin/bash

#$ -N fixVCF
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

raw_vcf=${1}
contigs="/hpcudd/home/boris/bin/dev/sofia/affy_scripts/contigs_hg19_for_vcf_header.txt"

awk '$0~/^##/' ${raw_vcf} > ${raw_vcf}.tmp.header.1
awk '$0~/^#CHROM/' ${raw_vcf} > ${raw_vcf}.tmp.header.2
awk '{gsub("\\.CEL","");print}' ${raw_vcf}.tmp.header.2 > ${raw_vcf}.tmp.header.3
awk '$0!~/^#/' ${raw_vcf} > ${raw_vcf}.unsorted.body

awk '$1~/[0-9]/' ${raw_vcf}.unsorted.body |sort -k1g,1g -k2g,2g > ${raw_vcf}.sorted.autosomes
awk '$1~/[XY]/' ${raw_vcf}.unsorted.body | sort -k1,1 -k2g,2g >  ${raw_vcf}.sorted.XY
awk '$1=="MT"' ${raw_vcf}.unsorted.body | sort -k2g,2g > ${raw_vcf}.sorted.MT

cat ${raw_vcf}.tmp.header.1 ${contigs} ${raw_vcf}.tmp.header.3 > ${raw_vcf}.header
cat ${raw_vcf}.header ${raw_vcf}.sorted.autosomes ${raw_vcf}.sorted.XY ${raw_vcf}.sorted.MT |bgzip -c > ${raw_vcf}.gz

tabix -p vcf ${raw_vcf}.gz

rm -f ${raw_vcf}.tmp.header.1 ${raw_vcf}.tmp.header.2 ${raw_vcf}.tmp.header.3 ${raw_vcf}.unsorted.body ${raw_vcf}.header ${raw_vcf}.sorted.autosomes ${raw_vcf}.sorted.XY ${raw_vcf}.sorted.MT

