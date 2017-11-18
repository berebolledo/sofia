#! /bin/bash

#$ -N gvcf-2
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

input=${1}
refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
bundle="/hpcudd/ICIM/shared/gatk-bundle"
dbsnp="${bundle}/b37/dbsnp_138.b37.vcf.gz"

gatk -Xms4g -Xmx8g \
    -T GenotypeGVCFs \
    -R $genome \
    -V $input \
    -o ${input}.raw.vcf

