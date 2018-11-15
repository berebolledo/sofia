#! /bin/bash

#$ -N metaCaller
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs


Gatk=${1}
Samtools=${2}
Freebayes=${3}
outname=${4}
model="/hpcudd/home/boris/bin/dev/sofia/vcf_prep/definitions.model"

/hpcudd/ICIM/boris/projects/gabriela/nwu/00-FINAL_VCF/all_raw_vcfs/raw-vcfs/VariantMetaCaller_v1.0/VariantMetaCaller \
    -c ${model} \
    -v HC HaplotypeCaller ${Gatk} \
    -v samtools samtools ${Samtools} \
    -v freebayes freebayes ${Freebayes}  |bgzip -c > ${outname}.metacaller.vcf.gz
