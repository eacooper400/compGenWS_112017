### Set the Library Path to Install and
### Load any Packages; install devtools
.libPaths(c('~/compGenWS_112017/Rlibs', .libPaths()))

### Set working directory
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/SampleData")
options(stringsAsFactors=FALSE)

### Read in the my functions file
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Read in the sample data
my.data=read.vcf("heli_phased.vcf", header=TRUE)

### Set the "anchor" SNP as 700,962
snp.of.interest = my.data[which(my.data$POS==700962),]

### Calculate R-squared between the SNP of
### interest and ALL other SNPs in the data
my.results = data.frame(Pos=my.data$POS, Rsq=0)
my.results$Rsq = apply(my.data, 1, function(x,y)
    calc_r2(as.vector(x, mode="character"), as.vector(y, mode="character")), y=snp.of.interest)

### Get the average LD in 10Kb bins
BinStarts=seq(1, max(my.results$Pos), 10000)
smoothed.results = data.frame(BinStarts, BinEnds=BinStarts+10000, AvgRsq=0)
for (i in 1:nrow(smoothed.results)) {
    w = subset(my.results, (my.results$Pos>=smoothed.results$BinStarts[i]
        & my.results$Pos<smoothed.results$BinEnds[i]))
    smoothed.results$AvgRsq[i] = mean(w$Rsq, na.rm=TRUE)
}

### Plot the results,
### add a star where the SNP of interest is
plot(smoothed.results$BinStarts, smoothed.results$AvgRsq,
     xlab="Position (bp)", ylab="Linkage Disequilibrium (r2)",
     type="l", col="red", ylim=c(0,1))
points(snp.of.interest[,2], 1, pch=8, col="red", cex=2)

### Now, pick a control SNP and do the same thing
control.snp = my.data[which(my.data$POS==1531867),]
my.results$Rsq2 = apply(my.data, 1, function(x,y)
    calc_r2(as.vector(x, mode="character"), as.vector(y, mode="character")), y=control.snp)
smoothed.results$AvgRsq2=0
for (i in 1:nrow(smoothed.results)) {
    w = subset(my.results, (my.results$Pos>=smoothed.results$BinStarts[i]
        & my.results$Pos<smoothed.results$BinEnds[i]))
    smoothed.results$AvgRsq2[i] = mean(w$Rsq2, na.rm=TRUE)
}
lines(smoothed.results$BinStarts, smoothed.results$AvgRsq2, col="blue")
points(control.snp[,2], 1, pch=8, col="blue", cex=1.5)
