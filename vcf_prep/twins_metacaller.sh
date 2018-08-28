#! /bin/bash

#$ -N metaCaller
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

chrom=${1}

Gatk="/hpcudd/ICIM/boris/projects/andresklein/twins/00_pipeline/gatk/bychr"
Samtools="/hpcudd/ICIM/boris/projects/andresklein/twins/00_pipeline/samtools/bychr"
Freebayes="/hpcudd/ICIM/boris/projects/andresklein/twins/00_pipeline/freebayes/bychr"
Metacaller="/hpcudd/ICIM/boris/projects/andresklein/twins/00_pipeline/metacaller"

/hpcudd/ICIM/boris/projects/gabriela/nwu/00-FINAL_VCF/all_raw_vcfs/raw-vcfs/VariantMetaCaller_v1.0/VariantMetaCaller \
    -c /hpcudd/ICIM/boris/projects/gabriela/nwu/00-FINAL_VCF/all_raw_vcfs/raw-vcfs/definitions.model \
    -v HC HaplotypeCaller ${Gatk}/gatk_twins_qual_exonic_chr_${chrom}.vcf \
    -v samtools samtools ${Samtools}/samtools_twins_qual_exonic_chr_${chrom}.vcf \
    -v freebayes freebayes ${Freebayes}/freebayes_twins_qual_exonic_chr_${chrom}.vcf > ${Metacaller}/metacaller_twins_qual_exonic_chr_${chrom}.vcf
