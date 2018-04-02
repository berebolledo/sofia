#! /bin/bash

#$ -N bird2vcf
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs


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

refdata="${genomes}/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"

specs="/hpcudd/ICIM/shared/GenomeWideSNP_6/specs/"
annotation="${specs}/GenomeWideSNP_6.na30.annot.hg19.csv.pickle.gz"


while getopts 'c:a:v:h' ARGS; do
	case "$ARGS" in
		c)
          calls="$OPTARG"
          ;;
        a)
          arrayid="$OPTARG"
          ;;
        v)
          vcfid="$OPTARG"
          ;;
        h)
          echo "script usage: $(basename $0) [-c calls.birdseed] [-a array.id] [-v vcf.id]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-c calls.birdseed] [-a array.id] [-v vcf.id]" >&2
          exit 1
          ;;
    esac
done


birdseed2vcf                             \
    --birdseed ${calls}                  \
    --output_vcf ${vcfid}.vcf            \
    --snp_annotation_file ${annotation}  \
    --array_sample ${arrayid}            \
    --vcf_sample ${vcfid}                \
    --fasta ${genome}                    \

