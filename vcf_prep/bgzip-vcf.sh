#! /bin/bash

#$ -N vcf.gz
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

bgzip ${1}
tabix -p vcf ${1}.gz
