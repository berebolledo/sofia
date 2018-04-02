#!/bin/bash

#$ -N tpm
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

file=${1}
sample=$(echo ${file}|cut -d "." -f 1)
reflen="/hpcudd/home/boris/projects/annmarie/definitive/gene_tx_len-gdc.txt"


join -j 1 \
    <(sort -k1,1 ${reflen}) \
    <(zcat ${file}|\
    awk '$1!~/__/ {gsub("\\.","\t");print $1"\t"$3}'|\
    sort -k1,1)|\
    awk 'BEGIN{OFS="\t"}{print $1,$2,($4/$3)*1000}' \
    > ${file}.tmp

sum=$(awk '{sum+=$3}END{print sum}' ${file}.tmp)
awk -v SUM=${sum} 'BEGIN{OFS="\t"}{print $1,$2,($3/SUM)*1000}' ${file}.tmp > ${sample}.tpm
rm -f ${file}.tmp
