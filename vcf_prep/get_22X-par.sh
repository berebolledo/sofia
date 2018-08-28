#! /bin/bash

#$ -N get22Xpar
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

input=${1}
bed="/hpcudd/home/boris/bin/dev/sofia/vcf_prep/ensembl_hg19_X-par.bed"

if [ ${input: -7} == ".vcf.gz" ]
then
    vcf=${input}
else
    bgzip ${input}
    vcf=${input}.gz
fi

vcftools              \
    --gzvcf ${vcf}    \
    --bed ${bed}      \
    --recode          \
    --recode-INFO-all \
    --stdout | bgzip -c > dip22Xpar_${1}

tabix -p vcf dip22Xpar_${1}
