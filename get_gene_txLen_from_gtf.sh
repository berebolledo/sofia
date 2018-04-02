awk -F"\t" '$3=="exon" {ID=substr($9, index($9,"transcript_id") + 15, 15); L[ID]+=$5-$4+1; GENE[ID]=substr($9, index($9,"gene_id") + 9, 15)} END{for(i in L){print GENE[i]"\t"i"\t"L[i]}}' genes.gtf
