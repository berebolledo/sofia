#! /bin/bash

#$ -N annovar
#$ -M brebolledo@udd.cl
#$ -m beas
#$ -o /hpcudd/home/boris/storage/data/logs
#$ -e /hpcudd/home/boris/storage/data/logs

input_vcf=${1}

zcat ${input_vcf}|awk 'BEGIN{OFS="\t"}{if($0~/^#/) print $0; else print $1,$2,$3,$4,$5,$6,$7,$8}' > ${input_vcf}-minimal.vcf

perl                                                                           \
    /hpcudd/home/boris/storage/data/variantAnnotation/annovar/table_annovar.pl \
    ${input_vcf}-minimal.vcf                                                   \
    /hpcudd/home/boris/storage/data/variantAnnotation/annovar/humandb          \
    -buildver hg19                                                             \
    -remove                                                                    \
    -out ${input_vcf}_annovar                                                  \
    -protocol refGene,esp6500siv2_all,1000g2015aug_all,avsnp144,dbnsfp30a,clinvar_20160302,exac03,dbscsnv11,dbnsfp31a_interpro,rmsk,ensGene,knownGene \
    -operation  g,f,f,f,f,f,f,f,f,r,g,g                                                                                                               \
    -nastring .                                                                                                                                       \
    -vcfinput

exitval=$?

if [ $exitval -eq 0 ]
then
  rm -f ${input_vcf}-minimal.vcf
fi

