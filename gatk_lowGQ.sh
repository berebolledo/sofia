#! /bin/bash

#$ -N lowGQ
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

input=${1}
refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"

gatk -Xms4g -Xmx8g \
    -T VariantFiltration \
    -R ${genome} \
    -V ${input} \
    -G_filter "GQ<20.0" \
    -G_filterName lowGQ \
    -o markedLowGQ_${input}
