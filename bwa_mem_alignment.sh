#! /bin/bash

#$ -N bwa
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe orte 8

read1=${1}
read2=${2}
RGID=${3}
RGLB=${4}

genomes="/hpcudd/ICIM/shared/genomes"
index="${genomes}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"

bwa mem                                                       \
    -t 8                                                      \
    -M                                                        \
    -R "@RG\tID:${RGID}\tLB:${RGLB}\tSM:${RGID}\tPL:ILLUMINA" \
    ${index}                                                  \
    ${read1}                                                  \
    ${read2}| samtools view -@ 2 -Sb -o ${RGID}.bam - 2>/dev/null 
