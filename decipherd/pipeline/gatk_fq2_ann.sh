#! /bin/bash

#$ -N fq2annvcf
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs
#$ -pe smp 8

#set -e
#set -u
#set -o pipefail


if [ $HOSTNAME == 'sofia.udd.cl' ] || [[ $HOSTNAME == compute*-1-*.local ]]
then
    genomes="/hpcudd/ICIM/shared/genomes"
    bundle="/hpcudd/ICIM/shared/gatk-bundle"
    annDir="/hpcudd/ICIM/boris/annotation"
    export PATH="/hpcudd/home/boris/miniconda3/bin:$PATH"
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

if [ $exit_bqsr1 -eq 0 ] && [ -s markDups.sorted.${RGID}_ALLchr_recal_data.table ]
then
    gatk -Xms4g -Xmx8g                                        \
        -T PrintReads                                         \
        -R $genome                                            \
        -I markDups.sorted.${RGID}.bam                        \
        -BQSR markDups.sorted.${RGID}_ALLchr_recal_data.table \
        -o bqsr.markDups.sorted.${RGID}.bam

    exit_bqsr2=$?
fi

if [ $exit_bqsr2 -eq 0 ] && [ -s bqsr.markDups.sorted.${RGID}.bam ]
then
    rm -f markDups.sorted.${RGID}.bam
    rm -f markDups.sorted.${RGID}.bai
    rm -f markDups.sorted.${RGID}_ALLchr_recal_data.table
    samtools view -b bqsr.markDups.sorted.${RGID}.bam|bamleftalign -f ${genome}|samtools view -b - > left-aligned.bqsr.markDups.sorted.${RGID}.bam
    samtools index left-aligned.bqsr.markDups.sorted.${RGID}.bam
    exit_left=$?
fi

if [ $exit_left -eq 0 ] && [ -s left-aligned.bqsr.markDups.sorted.${RGID}.bam ]
then
    rm -f bqsr.markDups.sorted.${RGID}.bam
    rm -f bqsr.markDups.sorted.${RGID}.bai

    gatk -Xms4g -Xmx8g     \
        -T HaplotypeCaller \
        -R $genome         \
        -I left-aligned.bqsr.markDups.sorted.${RGID}.bam \
        --genotyping_mode DISCOVERY      \
        --emitRefConfidence GVCF         \
        --variant_index_type LINEAR      \
        --variant_index_parameter 128000 \
        -o ${RGID}.gvcf
    exit_gvcf=$?
fi

if [ $exit_gvcf -eq 0 ] && [ -s ${RGID}.gvcf ]
then
    gatk -Xms4g -Xmx8g   \
        -T GenotypeGVCFs \
        -R $genome       \
        -V ${RGID}.gvcf  \
        -o ${RGID}.raw.vcf
    exit_vcf=$?
fi

if [ $exit_vcf -eq 0 ] && [ -s ${RGID}.raw.vcf ]
then
    python ${annDir}/InterVar/Intervar.py \
        --input=${RGID}.raw.vcf           \
        --input_type=VCF                  \
        --output=${RGID}.raw.vcf.intervar \
        --buildver=hg19                   \
        --database_intervar=${annDir}/InterVar/intervardb/           \
        --table_annovar=${annDir}/annovar/table_annovar.pl           \
        --annotate_variation=${annDir}/annovar/annotate_variation.pl \
        --database_locat=${annDir}/annovar/humandb/                  \
        --convert2annovar=${annDir}/annovar/convert2annovar.pl
    exit_intervar=$?
fi



if [ $exit_intervar -eq 0 ] && [ -s ${RGID}.raw.vcf.avinput ]
then

  input=${RGID}.raw.vcf
  python ../mergeAnnotations.py ${input}
  exit_end=$?

fi


if [ $exit_end -eq 0 ] && [ -s ${input}_annotations.tsv ]
then

    rm -f ${input}.avinput
    rm -f ${input}.${annExt}*
    rm -f ${input}.annovar.tmp.header
fi

echo -------------------------------
echo Exit codes
echo -------------------------------
echo bwa   = $exit_bwa
echo sort  = $exit_sort
echo dedup = $exit_mkd
echo bqsr1 = $exit_bqsr1
echo bqsr1 = $exit_bqsr2
echo left  = $exit_left
echo gvcf  = $exit_gvcf
echo vcf   = $exit_vcf
echo inter = $exit_intervar
echo final = $exit_end
echo -------------------------------