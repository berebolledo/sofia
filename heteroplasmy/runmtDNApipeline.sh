#! /bin/bash

file=${1}
ref='/storage/shared/references/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa'
readLength=`samtools view ${file} |head -100|awk '{sum+=length($10)} END{print sum/NR}'`

bamtools filter -region MT -isPaired true -isProperPair true -isMateMapped true -isPrimaryAlignment true -in ${file} -out proper.marked.srt.${file}

picard SortSam I=proper.marked.srt.${file} O=query.proper.marked.srt.${file} SO=queryname VALIDATION_STRINGENCY=SILENT

python rm_chim_in_pair.py query.proper.marked.srt.${file} ${readLength}

mv dechim.rlen.query.proper.marked.srt.${file} realigned.dechim.rlen.query.proper.marked.srt.${file}

picard SortSam I=realigned.dechim.rlen.query.proper.marked.srt.${file}  O=srt.realigned.dechim.rlen.query.proper.marked.srt.${file} SO=coordinate VALIDATION_STRINGENCY=SILENT

samtools index srt.realigned.dechim.rlen.query.proper.marked.srt.${file}

name2="srt.realigned.dechim.rlen.query.proper.marked.srt.${file}"

angsd/angsd -i ${name2} -doHetPlas 2 -out ${name2} -nThreads 4 -minQ 30 -minMapQ 20 -r MT -nLines 10000 -doCounts 1 -dumpCounts 3 -doQsDist 1 -howOften 10

python mitoMajorFromThorGL.py ${name2}.hetGL

samtools fillmd -b ${name2} ${name2}.hetGL.major.fa 2>/dev/null > md.${name2}

picard SortSam I=md.${name2}  O=query.md.${name2} SO=queryname VALIDATION_STRINGENCY=SILENT
python nm-ratio.select.py query.md.${name2}

name=$(echo ${file}|cut -d '.' -f 4)
mv nm-ratio.query.md.${name2} ${name}.nvcReady.bam

picard SortSam I=${name}.nvcReady.bam  O=srt.${name}.nvcReady.bam SO=coordinate VALIDATION_STRINGENCY=SILENT
samtools index srt.${name}.nvcReady.bam

python nvc/nvc/naive_variant_caller.py -b srt.${name}.nvcReady.bam -i srt.${name}.nvcReady.bam.bai -o ${name}.vcf -r ${ref}  -s -q 30 -m 20 -t uint32 --region MT &>/dev/null
python variant-annotator/allele-counts.py -i ${name}.vcf  -o ${name}.counts -n -s
