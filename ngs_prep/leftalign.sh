#! /bin/bash

#$ -N leftalign
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

ref="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"
samtools view -b ${1}|bamleftalign -f ${ref}|samtools view -b - > left-aligned.${1}
