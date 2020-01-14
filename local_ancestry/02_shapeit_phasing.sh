#! /bin/bash

ref="/home/boris/storage/00_papers/local_ancestry_pipeline/1000GP_Phase3"

/home/boris/storage/00_papers/local_ancestry_pipeline/shapeit/bin/shapeit \
    --input-ref ${ref}/1000GP_Phase3_chr${1}.hap.gz \
    ${ref}/1000GP_Phase3_chr${1}.legend.gz \
    ${ref}/1000GP_Phase3.sample \
    -B ${2}_chr${1} \
    --duohmm \
    --input-map ${ref}/genetic_map_chr${1}_combined_b37.txt \
    --exclude-snp ${2}_chr${1}.mendel.snp.strand.exclude \
    --output-max ${2}_chr${1}.haps.gz \
    ${2}_chr${1}.sample
