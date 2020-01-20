#!/usr/bin/env Rscript

setwd("/Users/boris/GoogleDrive/UDD/research/bioinformatics/SABIO/projects/andres_klein/anyelo")
library(qqman)
pvalues <- read.table("2020/pvalues_2020_new.txt",stringsAsFactors=F, header=T)
annotations <- read.table("snp_annotation.txt", stringsAsFactors=F, header=T)
data = cbind(annotations[,c("CHROMOSOME", "START", "SNP")], pvalues)
data$CHROMOSOME <- as.numeric(gsub('X', 20, data$CHROMOSOME))
rm(pvalues)

args <- commandArgs(trailingOnly=TRUE)
enz <- args[1]

#enz <- 'beta.D.Gal'
#enz <- 'Sphyngomyelin'
tmp <- data[data[[enz]]<=4.1e-06,c("SNP", enz)]
tmp$gene <-annotations[data[[enz]] <= 4.1e-06, "Symbol"]
tmp$chr <-annotations[data[[enz]] <= 4.1e-06, "CHROMOSOME"]
tmp$pos <-annotations[data[[enz]] <= 4.1e-06, "START"]
highlight <- tmp[!duplicated(tmp[c(2,3)]),]



png(paste("2020/", enz,'.plot.png',sep=''),width=6.5,height=6, units='in',res=300)
manhattan(data, chr="CHROMOSOME", bp="START", snp="SNP", p=enz, 
          chrlabs=c(1:19,'X'),  suggestiveline = -log10(4.1e-06),
          genomewideline = -log10(1.28e-08), ylim=c(0,10),
          col=c(rgb(26,26,26,100,maxColorValue=255),
                rgb(153,153,153,100,maxColorValue=255)),
          highlight=highlight$SNP, main=enz, cex.axis=0.4)
dev.off()

write.table(highlight, paste("2020/",enz,'.highlight.txt',sep=''), sep='\t',
            col.names=T, row.names=F, quote=F)
