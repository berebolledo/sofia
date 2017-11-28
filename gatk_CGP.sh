#! /bin/bash

#$ -N gcp
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

input=${1}

refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
bundle="/hpcudd/ICIM/shared/gatk-bundle/b37"
g1kp3="${bundle}/1000G_phase3_v4_20130502.sites.vcf.gz"

gatk -Xms4g -Xmx8g \
    -T CalculateGenotypePosteriors \
    -R ${genome} \
    -V ${input} \
    -supporting ${g1kp3} \
    -o posteriors_${input}
