#! /bin/bash

#$ -N varselection
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

input=${1}

if [ $2 == 'getindels' ]
then
    option='-selectType INDEL'
    vartype='indels'
elif [ $2 == 'rmindels' ]
then
    option='--selectTypeToExclude INDEL'
    vartype='not-indels'
else
    option='-selectType SNP'
    vartype='snps'
fi

refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"

gatk -Xms4g -Xmx8g \
    -T SelectVariants \
    -R $genome \
    -V $input \
    ${option} \
    -o ${vartype}.${input}

