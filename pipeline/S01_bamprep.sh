#! /bin/bash

#$ -N bamprep
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe orte 8

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

while getopts '1:2:i:l:h' ARGS; do
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
        h)
          echo "script usage: $(basename $0) [-1 read1.fq] [-2 read2.fq] [-i readgroup ID] [-l readgroup LIB]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-1 read1.fq] [-2 read2.fq] [-i readgroup ID] [-l readgroup LIB]" >&2
          exit 1
          ;;
    esac
done
    
shift "$(($OPTIND - 1))"
mkdir -p ${RGID}_tmpdir

index="${genomes}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"
refdata="${genomes}/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
dbsnp="${bundle}/b37/dbsnp_138.b37.vcf.gz"
indels="${bundle}/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"

bwa mem                                                       \
    -t 8                                                      \
    -M                                                        \
    -R "@RG\tID:${RGID}\tLB:${RGLB}\tSM:${RGID}\tPL:ILLUMINA" \
    ${index}                                                  \
    ${read1}                                                  \
    ${read2}| samtools view -@ 2 -Sb -o ${RGID}_tmpdir/${RGID}.bam - 2>/dev/null

exit_bwa=$?

cd ${RGID}_tmpdir 

if [ $exit_bwa -eq 0 ] && [ -s ${RGID}.bam ]
then
	picard SortSam                     \
        I=${RGID}.bam                \
        O=sorted.${RGID}.bam         \
        SO=coordinate                \
        VALIDATION_STRINGENCY=SILENT \
        TMP_DIR=tmpdir               \
        VERBOSITY=ERROR              \
        QUIET=true

    exit_sort=$?
fi


if [ $exit_sort -eq 0 ] && [ -s sorted.${RGID}.bam ]
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


if [ $exit_mkd -eq 0 ] && [ -s markDups.sorted.${RGID}.bam ]
then
    rm -f sorted.${RGID}.bam
    rm -fr tmpdir

    gatk -Xms4g -Xmx8g                 \
        -T BaseRecalibrator            \
        -R $genome                     \
        -I markDups.sorted.${RGID}.bam \
        -knownSites $dbsnp             \
        -knownSites $indels            \
        -o markDups.sorted.${RGID}_ALLchr_recal_data.table

    exit_bqsr1=$?
fi

if [ $exit_bqsr -eq 0 ] && [ -s markDups.sorted.${RGID}_ALLchr_recal_data.table ]
then
    gatk -Xms4g -Xmx8g                                        \
        -T PrintReads                                         \
        -R $genome                                            \
        -I markDups.sorted.${RGID}.bam                        \
        -BQSR markDups.sorted.${RGID}_ALLchr_recal_data.table \
        -o bqsr_markDups.sorted.${RGID}.bam

    exit_bqsr2=$?
fi

if [ $exit_bqsr1 -eq 0 ] && [ -s bqr_markDups.sorted.${RGID}.bam ]
then
    rm -f markDups.sorted.${RGID}.bam
    rm -f markDups.sorted.${RGID}.bai
fi