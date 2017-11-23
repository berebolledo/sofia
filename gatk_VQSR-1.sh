#! /bin/bash

#$ -N vqsr-1
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -cwd

source $HOME/.bashrc

input=${1}
mode=${2}

if [ ${mode} == "INDEL" ]
then
    extra_opt="--maxGaussians 4"
else
    extra_opt=""
fi

refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"

bundle="/hpcudd/ICIM/shared/gatk-bundle/b37"
hapmap="${bundle}/hapmap_3.3.b37.vcf.gz"
omni="${bundle}/1000G_omni2.5.b37.vcf.gz"
g1k="${bundle}/1000G_phase1.snps.high_confidence.b37.vcf.gz"
dbsnp="${bundle}/dbsnp_138.b37.vcf.gz"

gatk -Xms4g -Xmx8g \
    -T VariantRecalibrator \
    -R ${genome} \
    -input ${input} \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${g1k} \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode ${mode} \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile ${input}_recalibrate_${mode}.recal \
    -tranchesFile ${input}_recalibrate_${mode}.tranches \
    -rscriptFile ${input}_recalibrate_${mode}_plots.R \
    ${extra_opt}
