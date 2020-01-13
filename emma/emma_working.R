library(emma)

baseDir <- "/hpcudd/home/boris/storage/projects/andresklein"
setwd(paste0(baseDir, "/emma2020"))
load(paste0(baseDir, "/emma/2019/base_files/all_4million_98genotypes.RData"))

strains <- read.table("anyelo_strains.txt", header=F, stringsAsFactors=F)
colnames(strains) <- "strain"

data_strains <- allgtps[,c(colnames(allgtps)[1:4],strains$strain)]
rm(allgtps)

geno <- function(row){
  ref <- unlist(strsplit(row[4], split="/"))
  alleles <- row[5:31]
  transformed <- c()
  
  for (i in 1:length(alleles)){
    if (toupper(alleles[i]) == toupper(ref[1])){
      transformed <- c(transformed,0)}
    else if(toupper(alleles[i]) == toupper(ref[2])) {
      transformed <- c(transformed,1)}
    else {
      transformed <- c(transformed, NA)}
  }
  
  return(transformed)
}

geno_strains <- t(apply(data_strains,1,geno))
geno_strains <- data.frame(geno_strains)
colnames(geno_strains) <- strains$strain
rm(data_strains)

enzymes <- read.table("enzymes_corrected.txt", sep="\t", header=T, stringsAsFactors=F)
enzymes <- enzymes[enzymes$strain==strains$strain,]
ys <- as.matrix(enzymes[,2:17])
ys <- t(ys)
colnames(ys) <- strains$strain
xs <- as.matrix(geno_strains)
K <- emma.kinship(xs)

rs <- emma.ML.LRT(ys,xs,K)

pvalues <- rs$ps
colnames(pvalues) <- colnames(enzymes[2:17])

write.table(pvalues, "pvalues_new.txt",sep="\t",quote=F,row.names=F)
