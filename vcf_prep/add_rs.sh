#! /bin/bash

#$ -N add_rs
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

in_vcf=${1}
ref_vcf="/hpcudd/ICIM/shared/gatk-bundle/b37/dbsnp_138.b37.vcf.gz"

bcftools annotate            \
    --annotations ${ref_vcf} \
    --columns ID             \
    -O z -o rsid_${in_vcf}   \
    ${in_vcf}

tabix -p vcf rsid_${in_vcf}


