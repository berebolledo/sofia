#! /bin/bash

scripts="/home/boris/storage/00_papers/local_ancestry_pipeline/ancestry_pipeline"
ref="/home/boris/storage/00_papers/local_ancestry_pipeline/1000GP_Phase3"
example="/home/boris/storage/00_papers/local_ancestry_pipeline/ancestry_example"
#--ref_keep ${example}/ACB_example.ref
#--admixed_keep ${example}/ACB_example.notref

python ${scripts}/shapeit2rfmix.py \
    --shapeit_hap_ref ${2}_chr${1}.haps.gz \
    --shapeit_hap_admixed ${2}_chr${1}.haps.gz \
    --shapeit_sample_ref ${2}_chr${1}.sample \
    --shapeit_sample_admixed ${2}_chr${1}.sample \
    --ref_keep POP.ref \
    --admixed_keep CHL.notref \
    --chr ${1} \
    --genetic_map ${ref}/genetic_map_chr${1}_combined_b37.txt \
    --out CHL
