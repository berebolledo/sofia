#! /bin/bash

#$ -N stringtie
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe orte 8

bam=${1}
gtf="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/KRT5.gtf"

stringtie     \
    -G ${gtf} \
    -e -p 8   \
    -o ${bam}.quantification ${bam}
