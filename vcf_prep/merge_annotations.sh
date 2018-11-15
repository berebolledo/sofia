#! /bin/bash

#$ -N ANNOVARmerge
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

input_vcf=${1}
annotations=${2}

cat <(grep "#" ${annotations}) <(paste <(grep -v "#" ${annotations}) <(zcat ${input_vcf}|grep -v "#"|cut -f 9-))|bgzip -c > annovar_${input_vcf}
tabix -p vcf annovar_${input_vcf}
