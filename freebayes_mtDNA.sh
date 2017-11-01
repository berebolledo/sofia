#! /bin/bash

#$ -N freebayes
#$ -M brebolledo@udd.cl
#$ -m bes
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

refdata="/hpcudd/ICIM/shared/genomes/Homo_sapiens/Ensembl/GRCh37"
genome="${refdata}/Sequence/WholeGenomeFasta/genome.fa"
targets="exons_${2}.bed"

bamlist=${1}
output=${2}


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
    -t ${targets} |bgzip -c > ${output}.vcf.gz
