

### Read in the my functions file
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Then read in YOUR vcf file
### You will need to change the path and the name below!
my.data=read.vcf("~/genomics-workshop/SNPs/combined.vcf", header=TRUE, stringsAsFactors=FALSE)

### Convert data to a binary matrix
num.samples = ncol(my.data) - 9
new=matrix(0, nrow=nrow(my.data), ncol=num.samples)
colnames(new) = colnames(my.data)[10:(ncol(my.data))]
for (i in 1:nrow(my.data)) {
    r = as.vector(my.data[i,], mode="character")
    g = sapply(strsplit(r[10:length(r)], split=":"), `[[`, 1)
    g = gsub("\\/", "", g)
    g[g==".."] = NA
    g[g=="11"] = "2"
    g[g=="00"] = "0"
    g[g=="01"] = "1"
    new[i,] = as.numeric(g)
}
snp=t(new)
stree = nj(dist.gene(snp))
plot(stree)
