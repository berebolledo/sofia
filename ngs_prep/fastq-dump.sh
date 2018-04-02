#! /bin/bash

#$ -N fastq-dump
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

sra=${1}

fastq-dump --gzip --split-3 $sra
