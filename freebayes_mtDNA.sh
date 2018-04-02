#! /bin/bash

#$ -N freeb-MT
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

while getopts 'b:t:o:h' ARGS; do
    case "$ARGS" in
        b)
          bamlist="$OPTARG"
          ;;
        t)
          targets="$OPTARG"
          ;;
        o)
          output="$OPTARG"
          ;;
        h)
          echo "script usage: $(basename $0) [-b bam.list] [-t targets.bed] [-o output]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-b bam.list] [-t targets.bed] [-o output]" >&2
          exit 1
          ;;
    esac
done
    
shift "$(($OPTIND - 1))"

refdata="${genomes}/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"


#   -m --min-mapping-quality Q
#                   Exclude alignments from analysis if they have a mapping
#                   quality less than Q.  default: 1
#
#   -q --min-base-quality Q
#                   Exclude alleles from analysis if their supporting base
#                   quality is less than Q.  default: 0
#
#   -z --read-max-mismatch-fraction N
#                   Exclude reads with more than N [0,1] fraction of mismatches where
#                   each mismatch has base quality >= mismatch-base-quality-threshold
#                   default: 1.0
#
#   -C --min-alternate-count N
#                   Require at least this count of observations supporting
#                   an alternate allele within a single individual in order
#                   to evaluate the position.  default: 2
#
#   --min-coverage N
#                   Require at least this coverage to process a site. default: 0
#
#   -N --exclude-unobserved-genotypes
#                   Skip sample genotypings for which the sample has no supporting reads.
#
#   -G --min-alternate-total N
#                   Require at least this count of observations supporting
#                   an alternate allele within the total population in order
#                   to use the allele in analysis.  default: 1
#
#   -p --ploidy N   Sets the default ploidy for the analysis to N.  default: 2

freebayes             \
    -m 30             \
    -q 30             \
    -z 0.5            \
    -C 2              \
    --min-coverage 5  \
    -N                \
    --ploidy 1        \
    -G 2              \
    -L ${bamlist}     \
    -f ${genome}      \
    -t ${targets} > ${output}.vcf 
