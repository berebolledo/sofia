#! /bin/bash

#$ -N functVars
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

annotated_vcf=${1}

bcftools filter -i 'INFO/Func.refGene!~"^int" && INFO/Func.refGene!~"stream$"' ${annotated_vcf}|bgzip -c > ontarget_${annotated_vcf}
tabix -p vcf ontarget_${annotated_vcf}
