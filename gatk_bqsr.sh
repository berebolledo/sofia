#! /bin/bash

#$ -N bqsr
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

source $HOME/.bashrc

input=${1}
chrom=${2}
refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
bundle="/hpcudd/ICIM/shared/gatk-bundle"
dbsnp="${bundle}/b37/dbsnp_138.b37.vcf.gz"
indels="${bundle}/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"

gatk -Xms4g -Xmx8g \
    -T BaseRecalibrator \
    -R $genome \
    -I $input \
    -L $chrom \
    -knownSites $dbsnp \
    -knownSites $indels \
    -o ${input}_${chrom}_recal_data.table

