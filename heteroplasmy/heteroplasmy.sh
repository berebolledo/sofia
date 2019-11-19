#! /bin/bash

# From recentoy aligned BAM to MT nt counts

#scripts="/home/boris/storage/00_papers/paper_mtDNA/heteroplasmy/scripts"
scripts="/home/boris/bin/dev/sofia/heteroplasmy"
ref='/storage/shared/references/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
bam=${1}

# Average read lenght 
readLength=$(samtools view ${bam} |head -1000|awk '{sum+=length($10)} END{print sum/NR}')


# Get best alignments to MT
bamtools filter                    \
    -region MT                     \
    -isPaired true                 \
    -isProperPair true             \
    -isMateMapped true             \
    -isPrimaryAlignment true       \
    -in ${bam} -out proper.${bam}

# Sort by read name in preparation
picard SortSam                    \
    I=proper.${bam}               \
    O=query.proper.${bam}         \
    SO=queryname                  \
    VALIDATION_STRINGENCY=SILENT

# Remove chimeric reads
python ${scripts}/rm_chim_in_pair.py query.proper.${bam} ${readLength}
mv dechim.rlen.query.proper.${bam} checkpoint-1.${bam}

picard SortSam                    \
    I=checkpoint-1.${bam}         \
    O=srt.checkpoint-1.${bam}     \
    SO=coordinate                 \
    VALIDATION_STRINGENCY=SILENT

samtools index srt.checkpoint-1.${bam} 

prefix="srt.checkpoint-1.${bam}"

# Quick major/minor count for self-reference
angsd              \
    -i ${prefix}   \
    -doHetPlas 2   \
    -out ${prefix} \
    -nThreads 4    \
    -minQ 30       \
    -minMapQ 20    \
    -r MT          \
    -nLines 10000  \
    -doCounts 1    \
    -dumpCounts 3  \
    -doQsDist 1    \
    -howOften 10


# Generate self-reference
python ${scripts}/mitoMajorFromThorGL.py ${prefix}.hetGL

# Replace mismatch info based on self-reference
samtools \
    fillmd -b ${prefix} \
    ${prefix}.hetGL.major.fa 2>/dev/null > md.${prefix}

# Sort new bam
picard SortSam          \
    I=md.${prefix}      \
    O=query.md.${prefix} \
    SO=queryname        \
    VALIDATION_STRINGENCY=SILENT

# Remove reads with increased number of mismatches
python ${scripts}/nm-ratio.select.py query.md.${prefix}

# Rename final bam
name=$(echo ${bam}|cut -d "." -f 5)

mv nm-ratio.query.md.${prefix} ${name}.nvcReady.bam

# Sort bam and index
picard SortSam \
    I=${name}.nvcReady.bam  \
    O=srt.${name}.nvcReady.bam \
    SO=coordinate \
    VALIDATION_STRINGENCY=SILENT
samtools index srt.${name}.nvcReady.bam


# Run NVC
python ${scripts}/naive_variant_caller.py \
    -b srt.${name}.nvcReady.bam           \
    -i srt.${name}.nvcReady.bam.bai       \
    -o ${name}.vcf                        \
    -r ${ref}  -s -q 30 -m 20             \
    -t uint32 --region MT &>/dev/null

# Get the counts 
python ${scripts}/allele-counts.py \
    -i ${name}.vcf                 \
    -o ${name}.counts -n -s


# Tyding up
if [ -s ${name}.nvcReady.bam ] && [ -s ${name}.counts ]
then
	mkdir -p keep
	rm -f ${name}.nvcReady.bam
	mv ${bam} ${bam}.bai ${name}.counts ${name}.vcf keep/
	mv srt.${name}.nvcReady.bam srt.${name}.nvcReady.bam.bai keep/
	rm -f *.${bam}*
fi

