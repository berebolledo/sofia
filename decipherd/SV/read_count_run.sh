#! /bin/bash

#$ -N read_count
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

Rscript /hpcudd/home/boris/bin/dev/sofia/decipherd/read_count.R ${1} ${2} 
