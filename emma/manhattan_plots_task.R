#!/usr/bin/env Rscript

library(qqman)
library(tidyr)
library(calibrate)
################################################################################

myManhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
                                                                                  "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
                         genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
                         annotatePval = NULL, annotateTop = TRUE, ...) 
{
  CHR = BP = P = index = NULL
  if (!(chr %in% names(x))) 
    stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) 
    stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) 
    stop(paste("Column", p, "not found!"))
  if (!(snp %in% names(x))) 
    warning(paste("No SNP column found. OK unless you're trying to highlight."))
  if (!is.numeric(x[[chr]])) 
    stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
  if (!is.numeric(x[[bp]])) 
    stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) 
    stop(paste(p, "column should be numeric."))
  d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) 
    d = transform(d, SNP = x[[snp]])
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$logp <- -log10(d$P)
  }
  else {
    d$logp <- d$P
  }
  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  }
  else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index == 
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + 
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
                                                                     ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
                                            names(dotargs)]))
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  if (nchr == 1) {
    axis(1, ...)
  }
  else {
    axis(1, at = ticks, labels = labs, ...)
  }
  col = rep(col, max(d$CHR))
  if (nchr == 1) {
    with(d, points(pos, logp, pch = 20, col = col[1], ...))
  }
  else {
    icol = 1
    for (i in unique(d$index)) {
      with(d[d$index == unique(d$index)[i], ], points(pos, 
                                                      logp, col = col[icol], pch = 20, ...))
      icol = icol + 1
    }
  }
  if (suggestiveline) 
    abline(h = suggestiveline, col = "blue")
  if (genomewideline) 
    abline(h = genomewideline, col = "red")
  if (!is.null(highlight)) {
    if (any(!(highlight %in% d$SNP))) 
      warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col = "tomato", pch = 20, 
                             ...))
  }
  if (!is.null(annotatePval)) {
    topHits = subset(d, P <= annotatePval)
    par(xpd = TRUE)
    if (annotateTop == FALSE) {
      with(subset(d, P <= annotatePval), textxy(pos, -log10(P), 
                                                offset = 0.625, labs = topHits$SNP, cex = 0.45), 
           ...)
    }
    else {
      topHits <- topHits[order(topHits$P), ]
      topSNPs <- NULL
      for (i in unique(topHits$CHR)) {
        chrSNPs <- topHits[topHits$CHR == i, ]
        topSNPs <- rbind(topSNPs, chrSNPs[1, ])
      }
      textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
             labs = topSNPs$SNP, cex = 0.5, ...)
    }
  }
  par(xpd = FALSE)
}
################################################################################
myColor <- function(color){
  
  return(rgb(matrix(col2rgb(color),1,3), 
             alpha=10, maxColorValue=255))
}
################################################################################

setwd("/Users/boris/GoogleDrive/UDD/research/bioinformatics/SABIO/projects/andres_klein/anyelo")
pvalues <- read.table("2020/pvalues_2020_new.txt",stringsAsFactors=F, header=T)
annotations <- read.table("snp_annotation.txt", stringsAsFactors=F, header=T)
data <- cbind(annotations[,c("CHROMOSOME", "START", "Symbol")], pvalues)
data <-  cbind(data, unite(data[,1:3], "SNP", sep="_"))
data$CHROMOSOME <- as.numeric(gsub('X', 20, data$CHROMOSOME))
rm(pvalues)
gc()

args <- commandArgs(trailingOnly=TRUE)
enz <- args[1]

#enz <- 'beta.D.Gal'
enz <- 'GluCer'
tmp <- data[data[[enz]]<=4.1e-06,c("SNP", enz)]
tmp$gene <-annotations[data[[enz]] <= 4.1e-06, "Symbol"]
tmp$chr <-annotations[data[[enz]] <= 4.1e-06, "CHROMOSOME"]
tmp$pos <-annotations[data[[enz]] <= 4.1e-06, "START"]
highlight <- tmp[!duplicated(tmp[c(2,3)]),]

highlight <- highlight[order(highlight[[enz]]),]
highlight <- highlight[!duplicated(highlight$gene),]


yMax <- round(max(-log10(tmp[[enz]]))) + 1
col1 <- "steelblue4"
col2 <- "steelblue1"
tmp_data <- data[-log10(data[[enz]]) > 1 , c("CHROMOSOME", "START", "SNP", enz)]

png(paste("2020/", enz,'.plot.png',sep=''),width=6.5,height=6, units='in',res=300)
myManhattan(tmp_data, chr="CHROMOSOME", bp="START", snp="SNP", p=enz, 
          chrlabs=c(1:19,'X'),  suggestiveline = -log10(4.1e-06),
          genomewideline = -log10(1.28e-08), ylim=c(1,yMax),
          col=c(myColor(col1), myColor(col2)),
          highlight=highlight$SNP, main=enz, cex.axis=0.6)
dev.off()

write.table(highlight, paste("2020/",enz,'.highlight.txt',sep=''), sep='\t',
            col.names=T, row.names=F, quote=F)
