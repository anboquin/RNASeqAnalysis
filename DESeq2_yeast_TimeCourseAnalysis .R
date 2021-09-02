##Time course experiments - to find genes that react in a condition-specific manner over time 
#source: RNA-Seq workflow: gene-level exploratory analysis and differential expression

#here: demonstrate a basic time course analysis w the fission data package 
#package contains gene counts for an RNA-seq time course of fission yeast
#yeast were exposed to oxidative stress and half of the samples contain a deletion of the gene atf21 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install ("DESeq2")

library("DESeq2")
library("fission")

data("fission")
ddsTC <- DESeqDataSet(fission, ~strain + minute + strain:minute) 
      #formula models the strain difference at time 0, the difference over time, and any strain-specific differences over time

#perform a likelihood ratio test 
  #in this test two models are compared: 
      #a full model (includes all effects we are intersted in modelling)
      #a reduced model (removes specific effects we are interested in subjecting to null hypothesis testing)
  #here reduced model = remove interaction term bw strain and time (strain:minute)
ddsTC <- DESeq(ddsTC, test = "LRT", reduced = ~strain + minute)
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol
head(resTC[order(resTC$padj),],4)

#another option for time series data = model the counts as a smooth function of time & include an interaction term of the condition 
##plot count for the gene with the smallest adjusted pvalue
#testing for condition-dependent time profile and accounting for differences at time0 
data <- plotCounts(ddsTC, which.min(resTC$padj),
                   intgroup=c("minute","strain"), returnData=TRUE)
ggplot(data, aes(x=minute, y=count, color=strain, group=strain)) +
geom_point() + stat_smooth(se=FALSE,method="loess") + scale_y_log10()

#3investigate wald tests for log2fold changes at individual time points 
resultsNames(ddsTC) #get intercept 
res30 <- results(ddsTC, name="strainmut.minute30", test="Wald")
res30[which.min(resTC$padj),]

##cluster significant genes by their profiles 
#extract a matrix of the shrunken log2 fold changes 
betas <- coef(ddsTC)
colnames(betas)

#now plot the log2fold changes in a heatmap 
library("pheatmap")
topGenes <- head(order(resTC$padj),20)
mat <- betas[topGenes, -c(1,2)]
thr <- 3
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
         cluster_col=FALSE)


