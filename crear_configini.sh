vcf=${1}
out=${2}

echo "[InterVar]" >> config.ini.tmp
echo "buildver = hg19" >> config.ini.tmp 
echo "# hg19" >> config.ini.tmp 
echo "inputfile = ${vcf}" >> config.ini.tmp
echo "# the inputfile and the path  example/ex1.avinput hg19_clinvar_20151201.avinput" >> config.ini.tmp
echo "# tab-delimited will be better for including the other information" >> config.ini.tmp
echo "inputfile_type = VCF" >> config.ini.tmp
echo "# the input file type VCF(vcf file with single sample),AVinput,VCF_m(vcf file with multiple samples)" >> config.ini.tmp
echo "outfile = ${out}" >> config.ini.tmp

cat config.ini.tmp parte_fija.tmp > config.ini
