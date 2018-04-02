#! /bin/bash

#$ -N cvrg_angsd
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

set -e
set -u
set -o pipefail

while getopts 'b:t:o:h' ARGS; do
        case "$ARGS" in
        b)
          bam="$OPTARG"
          ;;
        t)
          targets="$OPTARG"
          ;;
        o)
          out="$OPTARG"
          ;;
        h)
          echo "script usage: $(basename $0) [-b bam] [-t targets.bed] [-o prefix]" >&2
          exit 0
          ;;
        ?)
          echo "script usage: $(basename $0) [-b bam] [-t targets.bed] [-o prefix]" >&2
          exit 1
          ;;
    esac
done
    
shift "$(($OPTIND - 1))"




angsd               \
    -doCounts 1     \
    -dumpCounts 1   \
    -rf ${targets}  \
    -i ${bam}       \
    -out ${out}


nreads=$(samtools view -c ${bam})
cvrg_1=$(zcat ${out}.pos.gz|awk 'NR>1 && $3>0 {sum+=1} END{print sum}')
cvrg_3=$(zcat ${out}.pos.gz|awk 'NR>1 && $3>=3 {sum+=1} END{print sum}')
cvrg_5=$(zcat ${out}.pos.gz|awk 'NR>1 && $3>=5 {sum+=1} END{print sum}')
not_cvrg=$(expr 16569 - $cvrg_1)

mean_cvrg=$(zcat ${out}.pos.gz|awk 'NR>1 {sum+=$3} END{print sum/16569}')
median_cvrg=$(zcat ${out}.pos.gz|tail -n +2|cut -f 3|cat - <( printf '0%.0s\n' $(seq 1 $not_cvrg) )|datamash median 1)

over_mean=$(zcat ${out}.pos.gz|awk -v mean=${mean_cvrg} 'NR>1 && $3>=mean {sum+=1} END{print sum}')
over_median=$(zcat ${out}.pos.gz|awk -v median=${median_cvrg} 'NR>1 && $3>=median {sum+=1} END{print sum}')

rm -f ${out}.metrics
touch ${out}.metrics
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" sample nreads cvrg1 cvrg3 cvrg5 mean_cvrg median_cvrg bp_over_mean bp_over_median >> ${out}.metrics
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $out $nreads $cvrg_1 $cvrg_3 $cvrg_5 $mean_cvrg $median_cvrg $over_mean $over_median >> ${out}.metrics



