#! /bin/bash

#$ -N apt-geno
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

set -e
set -u
set -o pipefail

while getopts 'c:o:h' ARGS; do
	case "$ARGS" in
        c)
          cel_files="$OPTARG"
          ;;
        o)
          output="$OPTARG"
          ;;
        h)
          echo "script usage: $(basename $0) [-c celfiles.list] [-o output.prefix]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-c celfiles.list] [-o output.prefix]" >&2
          exit 1
          ;;
    esac
done
    
shift "$(($OPTIND - 1))"



first_line=$(head -1 ${cel_files})
lib="/hpcudd/ICIM/shared/GenomeWideSNP_6/libfiles"

if [ ${first_line} != "cel_files" ]
then
    cat <(echo "cel_files") ${cel_files} > ${cel_files}.tmp
fi

apt-probeset-genotype                          \
    -o ${output}                               \
    -a birdseed-dev                            \
    -c ${lib}/GenomeWideSNP_6.cdf                     \
    --set-gender-method cn-probe-chrXY-ratio          \
    --chrX-probes ${lib}/GenomeWideSNP_6.chrXprobes   \
    --chrY-probes ${lib}/GenomeWideSNP_6.chrYprobes   \
    --special-snps ${lib}/GenomeWideSNP_6.specialSNPs \
    --read-models-birdseed ${lib}/GenomeWideSNP_6.birdseed-v2.models \
    --cel-files ${cel_files}

if [ -s ${cel_files}.tmp ]
then
    rm -f ${cel_files}.tmp
fi