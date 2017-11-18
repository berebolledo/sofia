#! /bin/bash

#$ -N gvcf
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
    -T HaplotypeCaller \
    -R $genome \
    -I $input \
    --genotyping_mode DISCOVERY \
    --emitRefConfidence GVCF \
    -o ${input}.g.vcf

