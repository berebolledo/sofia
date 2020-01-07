library(Rsubread)

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
prefix <- args[2]
outputFile <- paste0(prefix, "_counts.txt")
bedFile <- "./refGene.hg19.sorted.nochr.bed"
bed <- read.table(bedFile, sep="\t", header=F, stringsAsFactors=F)
df <- data.frame(cbind(1:nrow(bed), bed[,c(1:3,6)]))
colnames(df) <- c("GeneID","Chr", "Start", "End","Strand")

res <-  featureCounts(files=file, annot.ext=df, allowMultiOverlap=TRUE, isPairedEnd=TRUE)


counts <- data.frame(res$counts)
colnames(counts) <- "count"
annotation <- data.frame(res$annotation)
out <- cbind(annotation, counts) 

total_reads <- sum(res$stat[2])

pm <- total_reads / 1000000
rpm <-  as.numeric(out$count) / pm
out$fpkm <- rpm / (as.numeric(out$Length) / 1000)


rpk <- out$count / (as.numeric(out$Length) / 1000)
pm <- sum(rpk) / 1000000
out$tpm <- rpk/pm

write.table(out, file=outputFile, sep="\t", col.names=T, row.names=F, quote=F)