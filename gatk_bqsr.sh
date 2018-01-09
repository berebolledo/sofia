#! /bin/bash

#$ -N bqsr
#$ -M brebolledo@udd.cl
#$ -m bes
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

while getopts 'i:c:h' ARGS; do
	case "$ARGS" in
        i)
          input="$OPTARG"
          ;;
        c)
          chrom="$OPTARG"
          ;;
        h)
          echo "script usage: $(basename $0) [-i input.bam] [-c chromosome N|ALL]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-i input.bam] [-c chromosome N|ALL]" >&2
          exit 1
          ;;
    esac
done
    
shift "$(($OPTIND - 1))"

refdata="${genomes}/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
dbsnp="${bundle}/b37/dbsnp_138.b37.vcf.gz"
indels="${bundle}/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"

if [ $chrom == "ALL" ]
then
    gatk -Xms4g -Xmx8g \
        -T BaseRecalibrator \
        -R $genome \
        -I $input \
        -knownSites $dbsnp \
        -knownSites $indels \
        -o ${input}_${chrom}_recal_data.table
else
    gatk -Xms4g -Xmx8g \
        -T BaseRecalibrator \
        -R $genome \
        -I $input \
        -L $chrom \
        -knownSites $dbsnp \
        -knownSites $indels \
        -o ${input}_${chrom}_recal_data.table
fi
