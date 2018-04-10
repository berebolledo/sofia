#! /bin/bash

#$ -N fixbirdseed
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

header="/hpcudd/home/boris/bin/dev/sofia/vcf_prep/header_v42_hg19_birdseed.txt"


while getopts 'v:h' ARGS; do
	case "$ARGS" in
        v)
          vcf_raw="$OPTARG"
          ;;
        h)
          echo "script usage: $(basename $0) [-v raw.vcf]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-v raw.vcf]" >&2
          exit 1
          ;;
    esac
done

head -4 ${vcf_raw}|tail -3 > ${vcf_raw}.header
awk '$1~/^[1-9]/' ${vcf_raw} > ${vcf_raw}.autosomes
awk '$1=="X"' ${vcf_raw} > ${vcf_raw}.X
awk '$1=="Y"' ${vcf_raw} > ${vcf_raw}.Y
awk '$1=="MT"' ${vcf_raw} > ${vcf_raw}.MT

sort -k1g,1g -k2g,2g ${vcf_raw}.autosomes > ${vcf_raw}.autosomes.sorted
sort -k2g,2g ${vcf_raw}.X > ${vcf_raw}.X.sorted
sort -k2g,2g ${vcf_raw}.Y > ${vcf_raw}.Y.sorted
sort -k2g,2g ${vcf_raw}.MT > ${vcf_raw}.MT.sorted

cat ${header} ${vcf_raw}.header ${vcf_raw}.autosomes.sorted ${vcf_raw}.X.sorted ${vcf_raw}.Y.sorted ${vcf_raw}.MT.sorted > fixed.${vcf_raw}

bgzip fixed.${vcf_raw}
tabix -p vcf fixed.${vcf_raw}.gz

rm -f ${vcf_raw}.header ${vcf_raw}.autosomes* ${vcf_raw}.X* ${vcf_raw}.Y* ${vcf_raw}.MT*
