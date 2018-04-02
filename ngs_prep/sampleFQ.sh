#! /bin/bash

#$ -N fqsample
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

fastq=${1}
fraction=${2}

seqtk sample ${fastq} ${fraction}|gzip -c > sample_${fraction}_${fastq}

