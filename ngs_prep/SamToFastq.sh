#! /bin/bash

#$ -N Bam2Fq
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

inputbam=${1}
outputname=${inputbam%.*}

picard SamToFastq                   \
    I=${inputbam}                   \
    F=${outputname}_1.fq.gz         \
    F2=${outputname}_2.fq.gz        \
    FU=${outputname}_unpaired.fq.gz \
    VALIDATION_STRINGENCY=SILENT    \
    TMP_DIR=${outputname}_tmpdir    \
    VERBOSITY=ERROR                 \
    QUIET=true

# Clean up

rm -fr ${outputname}_tmpdir

if [ ! -s ${outputname}_unpaired.fq.gz ] ; then
    rm -f ${outputname}_unpaired.fq.gz
fi
