#! /bin/bash

cd /home/boris/storage/00_papers/local_ancestry_pipeline/RFMix_v1.5.4

python RunRFMix.py \
    -e 2 -w 0.2 \
    --num-threads 2 \
    --use-reference-panels-in-EM \
    --forward-backward PopPhased \
    mydata/${2}_chr${1}.alleles \
    mydata/${2}.classes \
    mydata/${2}_chr${1}.snp_locations \
    -o mydata/${2}_chr${1}.rfmix
