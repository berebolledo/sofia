#! /bin/bash

#$ -N INFO_vcf
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

vcf=${1}
info=${2}


if [[ "$vcf" == *.gz$ ]]
then
    vcfopt="--gzvcf"
else
    vcfopt="--vcf"
fi

vcftools \
    ${vcfopt} ${vcf} \
    --get-INFO ${info} \
    --out ${vcf}_${info}.txt
