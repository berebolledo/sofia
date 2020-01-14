#! /bin/bash

scripts="/home/boris/storage/00_papers/local_ancestry_pipeline/ancestry_pipeline"

python ${scripts}/collapse_ancestry.py \
    --rfmix ${2}_chr${1}.rfmix.2.Viterbi.txt \
    --snp_locations ${2}_chr${1}.snp_locations \
    --fbk ${2}_chr${1}.rfmix.2.ForwardBackward.txt \
    --fbk_threshold 0.9 \
    --ind ${3} \
    --ind_info CHL.sample \
    --pop_labels EUR,AFR,NAT \
    --out ${3}

