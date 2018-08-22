#! /bin/bash

#$ -N BamIndex
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs


inputbam=${1}

picard BuildBamIndex             \
    I=${inputbam}                \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=${inputbam}_tmpdir \
    VERBOSITY=ERROR              \
    QUIET=true

# Clean up
rm -fr ${inputbam}_tmpdir
