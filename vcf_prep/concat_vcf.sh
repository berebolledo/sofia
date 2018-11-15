#! /bin/bash

#$ -N bcfconcat
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

file_list=${1}
output=${2}

bcftools concat     \
    -a -D -O z      \
    -f ${file_list} \
    -o ${output}.vcf.gz

tabix -p vcf ${output}.vcf.gz
