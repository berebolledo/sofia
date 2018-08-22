#! /bin/bash

#$ -N sam_vcf
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs

ref="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"
bamlist=${1}
pheno=${2}
out=${3}

bcftools mpileup -q 30 -Q 30 -d 1000000 -Ou -f ${ref} -b ${bamlist} | bcftools call --ploidy GRCh37 -S ${pheno} -vmO z -o ${out}.vcf.gz
