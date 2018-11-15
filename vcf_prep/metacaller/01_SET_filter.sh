#! /bin/bash

#$ -N filter_1
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

# SET filter: at least one allele was found by all three callers

bcftools filter -i 'INFO/SET~"HC+samtools+freebayes"' -O z -o SET-filter_${1} ${1}
tabix -p vcf SET-filter_${1}
