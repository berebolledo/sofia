#! /bin/bash

#$ -N wget
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

url=${1}
out="${url##*/}"

wget -q -O ${out} ${url} 
