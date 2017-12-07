#! /bin/bash

#$ -N bwa
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe orte 8
#$ -cwd

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


if [ $HOSTNAME == 'sofia.udd.cl' ]; then
    genomes="/hpcudd/ICIM/shared/genomes"
elif [ $HOSTNAME == 'mendel' ]; then
    genomes="/storage/shared/references"
else
    echo "Unrecognized host $HOSTNAME"
    echo "can't locate genome references"
    exit 1
 fi


mkdir -p ${RGID}_tmpdir
index="${genomes}/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa"


bwa mem                                                       \
    -t 8                                                      \
    -M                                                        \
    -R "@RG\tID:${RGID}\tLB:${RGLB}\tSM:${RGID}\tPL:ILLUMINA" \
    ${index}                                                  \
    ${read1}                                                  \
    ${read2}| samtools view -@ 2 -Sb -o ${RGID}_tmpdir/${RGID}.bam - 2>/dev/null