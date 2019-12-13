#! /bin/bash

input=${1}
name=$(echo $input|cut -d "." -f 1)

perl /storage/annovar/table_annovar.pl \
     ${input} \
     /storage/annovar/humandb  \
    -buildver hg19 -remove     \
    -out ${input}.annovar.tmp  \
    -protocol refGene,1000g2015aug_all,exac03,esp6500siv2_all,avsnp147,clinvar_20190305,dbnsfp33a,gnomad_exome,gnomad_genome  \
    -operation g,f,f,f,f,f,f,f,f -nastring . -vcfinput

n=$(head -1 ${input}.annovar.tmp.hg19_multianno.txt|datamash transpose|wc -l)
t=$(tail -n +2 ${input}.annovar.tmp.hg19_multianno.txt|head -1|datamash transpose|wc -l)
x=$(expr ${t} - ${n})
paste <(head -1 ${input}.annovar.tmp.hg19_multianno.txt) <(for i in $(seq 1 $x);do echo "Otherinfo";done|datamash transpose) > ${input}.annovar.tmp.header
cat ${input}.annovar.tmp.header <(tail -n +2 ${input}.annovar.tmp.hg19_multianno.txt) > ${input}_annovar.txt
cp ${input}_annovar.txt ${name}_annovar.txt

rm -f ${input}.annovar.tmp.*


python /storage/InterVar/Intervar.py \
--input=${input} \
--input_type=VCF \
--output=${input}.intervar \
--buildver=hg19 \
--database_intervar=/storage/InterVar/intervardb \
--table_annovar=/storage/annovar/table_annovar.pl \
--convert2annovar=/storage/annovar/convert2annovar.pl \
--annotate_variation=/storage/annovar/annotate_variation.pl \
--database_locat=/storage/annovar/humandb

mv ${input}.intervar.hg19_multianno.txt.intervar ${input}_intervar.txt
cp ${input}_intervar.txt ${name}_intervar.txt

rm -f ${input}.avinput
rm -f ${input}.intervar.hg19_multianno.txt
rm -f ${input}.intervar.hg19_multianno.txt.grl_p
