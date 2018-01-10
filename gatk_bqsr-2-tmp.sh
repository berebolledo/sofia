#! /bin/bash

#$ -N bqsr-2
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

input=${1}
chrom=${2}
refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"

gatk -Xms4g -Xmx8g \
    -T PrintReads \
    -R $genome \
    -I $input \
    -BQSR ${input}_${chrom}_recal_data.table \
    -o ${input}_bqsr.bam

