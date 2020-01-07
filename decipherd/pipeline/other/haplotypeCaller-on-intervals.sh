#! /bin/bash

#$ -N haploCall
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

source $HOME/.bashrc

set -e
set -u
set -o pipefail


if [ $HOSTNAME == 'sofia.udd.cl' ] || [[ $HOSTNAME == compute-1-*.local ]]
then
    genomes="/hpcudd/ICIM/shared/genomes"
    bundle="/hpcudd/ICIM/shared/gatk-bundle"
elif [ $HOSTNAME == 'mendel' ]
then
    genomes="/storage/shared/references"
    bundle="/storage/shared/gatk-bundle"
else
    echo "Unrecognized host $HOSTNAME"
    echo "can't locate genome references"
    exit 1
fi

while getopts 'i:t:h' ARGS; do
	case "$ARGS" in
        i)
          input="$OPTARG"
          ;;
        t)
          target="$OPTARG"
          ;;  
        h)
          echo "script usage: $(basename $0) [-i input.bam] [-t target]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-i input.bam] [-t target]" >&2
          exit 1
          ;;
    esac
done

refdata="${genomes}/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
intervals="${refdata}/Annotation/Intervals"

padding="--interval_padding 100"

if [ $target == 'auto-exons' ]
then
    interval=${intervals}/autosomes_exons.bed
    
elif [ $target == 'auto-whole' ]
then
    interval=${intervals}/autosomes.list

elif [ $target == 'X-exons' ]
then
    interval=${intervals}/X_exons.bed

elif [ $target == 'X-whole' ]
then
    interval="X"

elif [ $target == 'male-exons' ] 
then
    interval=${intervals}/Y_exons.bed

elif [ $target == 'male-XY' ]
then
    interval="-L X -L Y"


    
gatk -Xms4g -Xmx8g \
    -T HaplotypeCaller \
    -R $genome \
    -I $input \
    --genotyping_mode DISCOVERY \
    --emitRefConfidence GVCF \
    -o ${input}.g.vcf
