#!/bin/bash

vcf=${1}

cat ${vcf}|bcftools filter -e 'INFO/Func.refGene ~ "int" || INFO/Func.refGene ~ "down" || INFO/Func.refGene ~ "up"' > exonic_${vcf}
