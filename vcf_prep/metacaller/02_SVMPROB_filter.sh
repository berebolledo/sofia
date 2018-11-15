#! /bin/bash

#$ -N filter_2
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

# SVMPROB filter
# sum(SVMPROB) >= 0.5 (at least one allele >= 0.5) 

bcftools filter -i 'SUM(INFO/SVMPROB)>=0.5' -O z -o SVMPROB-filter_${1} ${1}
tabix -p vcf SVMPROB-filter_${1}
