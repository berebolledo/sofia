#! /bin/bash

#$ -N vqsr-2
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

input=${1}
mode=${2}

refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"

gatk -Xms4g -Xmx8g \
    -T ApplyRecalibration \
    -R ${genome} \
    -input ${input} \
    -mode ${mode} \
    --ts_filter_level 99.0 \
    -recalFile ${input}_recalibrate_${mode}.recal \
    -tranchesFile ${input}_recalibrate_${mode}.tranches \
    -o recalibrated_${mode}_${input}
