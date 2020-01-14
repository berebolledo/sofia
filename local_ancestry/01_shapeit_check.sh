#! /bin/bash

ref="/home/boris/storage/00_papers/local_ancestry_pipeline/1000GP_Phase3"
/home/boris/storage/00_papers/local_ancestry_pipeline/shapeit/bin/shapeit \
    -check \
    --input-ref ${ref}/1000GP_Phase3_chr${1}.hap.gz \
    ${ref}/1000GP_Phase3_chr${1}.legend.gz \
    ${ref}/1000GP_Phase3.sample \
    -B ${2}_chr${1} \
    --input-map ${ref}/genetic_map_chr${1}_combined_b37.txt \
    --output-log ${2}_chr${1}.mendel

