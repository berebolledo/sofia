#! /bin/bash

#$ -N vcf2plink
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

vcf=${1}
prefix=${2}

plink2           \
    --double-id  \
    --vcf ${vcf} \
    --make-bed   \
    --out ${prefix}
