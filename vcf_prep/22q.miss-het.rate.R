data = read.table("22q.miss-het.rate.txt",header=F, stringsAsFactors=F)
colnames(data) = c("sample", "heterozygosity", "missingness")
plot(
  log10(data$missingness), 
  data$heterozygosity,
  ylab="sample heterozygosityrate",
  xlab="sample missingness",
  xlim=c(-3,-2), 
  ylim=c(0.1,0.2),
  axes=F,
  bty="n",
  pch=21,
  col="blue",
  bg=rgb(100,149,237,150,maxColorValue=255))

axis(1, at=c(-3,-2), lab=c("0.1%","1%"))
axis(2, at=c(0.1, 0.15, 0.2), las=2)

meanHet  = mean(data$heterozygosity)
sdHet    = sd(data$heterozygosity)
upperlim = meanHet+3*sdHet
lowerlim = meanHet-3*sdHet

abline(h=upperlim, lty=2, col="tomato")
abline(h=lowerlim, lty=2, col="tomato")

for (i in 1:nrow(data)){
  if (data[i,2]>=upperlim | data[i,2]<=lowerlim){
    text(log10(data[i,3]),data[i,2],labels=data[i,1], pos=4, cex=0.75)
  }
}
