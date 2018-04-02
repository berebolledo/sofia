#! /bin/bash

#$ -N bamFilter
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

bamtools filter -in ${1} -out ${2}_${1} -region ${2}
