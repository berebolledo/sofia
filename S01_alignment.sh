#! /bin/bash

#$ -N bwa
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe orte 8

set -e
set -u
set -o pipefail

genomes="/hpcudd/ICIM/shared/genomes"
index="${genomes}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"

while getopts '1:2:i:l:' ARGS; do
	case "$ARGS" in
		1)
          read1="$OPTARG"
          ;;
        2)
          read2="$OPTARG"
          ;;
        i)
          RGID="$OPTARG"
          ;;
        l)
          RGLB="$OPTARG"
          ;;
        ?)
          echo "script usage: $(basename $0) [-1 read1.fq] [-2 read2.fq] [-i readgroup ID] [-l readgroup LIB]" >&2
          exit 1
          ;;
    esac
done    
shift "$(($OPTIND - 1))"

bwa mem                                                       \
    -t 8                                                      \
    -M                                                        \
    -R "@RG\tID:${RGID}\tLB:${RGLB}\tSM:${RGID}\tPL:ILLUMINA" \
    ${index}                                                  \
    ${read1}                                                  \
    ${read2}| samtools view -@ 2 -Sb -o ${RGID}_tmpdir/${RGID}.bam - 2>/dev/null

exit_bwa=$?

mkdir -p ${RGID}_tmpdir && cd ${RGID}_tmpdir 

if [ $exit_bwa -eq 0 && -s ${RGID}.bam  ]
then
	picard SortSam                   \
        I=${RGID}.bam                \
        O=sorted.${RGID}.bam         \
        SO=coordinate                \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=tmpdir               \
        VERBOSITY=ERROR              \
        QUIET=true

    exit_sort=$?
fi


if [ $exit_sort -eq 0 && -s sorted.${RGID}.bam ]
then
    rm -f ${RGID}.bam

    picard MarkDuplicates              \
        I=sorted.${RGID}.bam           \
        O=markDups.sorted.${RGID}.bam  \
        M=${RGID}.metrics.txt          \
        ASO=coordinate                 \
        VALIDATION_STRINGENCY=SILENT   \
        TMP_DIR=tmpdir                 \
        VERBOSITY=ERROR                \
        QUIET=true                     \
        CREATE_INDEX=true

    exit_mkd=$?
fi


if [ $exit_mkd -eq 0 && -s markDups.sorted.${RGID}.bam ]
then
    rm -f sorted.${RGID}.bam
    rm -fr tempdir
fi