#from: ;RNA-Seq workflow: gene-level exploratory analysis and differential expression (Love et al, 2015)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install ("DESeq2")
library("airway")
library("DESeq2")

#extdata is where R packages store external data 
dir <- system.file("extdata", package = "airway", mustWork = TRUE)
list.files(dir) #eight BAM files 

#load files
csvfile <- file.path(dir, "sample_table.csv")
(sampleTable <- read.csv(csvfile, row.names = 1)) #() trick to print the output 

##generate count matrix using summarizeOverlaps method 
#construct full paths to the files we want to perform counting operation on
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames) 

library("Rsamtools") #to provide an R interface to BAM files 
bamfiles <- BamFileList(filenames, yieldSize = 2000000) #specify how BAM files should be treated 

seqinfo(bamfiles[1]) #check chromosome names 

#defining gene models 
library("GenomicFeatures")
gtffile <- file.path(dir, "Homo_sapiens.GRCh37.75_subset.gtf")
(txdb <- makeTxDbFromGFF(gtffile, format="gtf", circ_seqs = character())) #txdb = transcript database

#produce GRangesList of all exons grouped by gene 
(ebg <- exonsBy(txdb, by="gene"))


library("GenomicAlignments")
#library("BiocParallel") if using multiple cores, this can sspeed up the counting process 
#register(SerialParam())

#create summarizedexperiment with counts 
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode = "Union", #describes what kind of read overlaps will be counted
                        singleEnd = FALSE,
                        ignore.strand = TRUE,
                        fragments = TRUE)
se

rowRanges(se)
str(metadata(rowRanges(se))) #display metadata compactly 

#assign the sampleTable as the colData of the se by converting it to a dataframe
(colData(se) <- DataFrame(sampleTable))


#Starting from SummarizedExperiment
data("airway")
se <- airway

#specify that untrt is the reference level for the dex variable 
se$dex <- relevel(se$dex, "untrt")
se$dex

#check the millions of fragments that uniquely aligned to the genes (2nd argument specifies decimal points)
round(colSums(assay(se))/1e6, 1)

colData(se)

dds <- DESeqDataSet(se, design = ~cell + dex)  #test for the effect of dexamethasone controlling for the effect of different donors’ cells
dds


#starting from count matrices
countdata <- assay(se)
head(countdata, 3)

coldata <- colData(se)

#construct deseqdataset object
(ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = coldata,
                                  design = ~cell + dex))

## exploratory analysis and visualization 
#prefiltering
nrow(dds) #64102
dds <- dds[rowSums(counts(dds)) >1, ]
nrow(dds) #29391

#rlog transformation 
rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

#to show effect of transformatiion plot 1st sample against 2n using the log2 function
par( mfrow = c( 1, 2 ) )
dds <- estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1),
     pch=16, cex=0.3)
plot(assay(rld)[,1:2],
     pch=16, cex=0.3)
  #we can see the rlog transformation compresses differences for the low count genes for which the data provide little info about differential expression

#sample distances - to assess overall similarity bw samples
sampleDists <- dist(t(assay(rld))) #t = transpose 
sampleDists

#visualize distance in a heatmap
library("pheatmap")
library("RColorBrewer")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$dex, rld$cell, sep="-") #changed rownames of the dist matrix to contain treatment type and patient number  
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors)

#another way to calculaate sample distances = Poisson distance 
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

#plot poisson distance heatmap
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(rld$dex, rld$cell, sep="-")
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

#pca plot 
plotPCA(rld, intgroup = c("dex", "cell"))

#building pca plot from scratch using ggplot2
#ask plotPCA function to return the data used for plotting
(pcadata <- plotPCA(rld, intgroup = c( "dex", "cell"), returnData=TRUE))
percentVar <- round(100 * attr(pcadata, "percentVar"))

library("ggplot2")
ggplot(pcadata, aes(PC1, PC2, color = dex, shape = cell)) + 
  geom_point(size=3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

#MDS plot - useful when we dont have a matrix of data but only a matrix of distances 
mdsdata <- data.frame(cmdscale(sampleDistMatrix))
mds <- cbind(mdsdata, as.data.frame(colData(rld)))
ggplot(mds, aes(X1, X2, color=dex, shape=cell)) + geom_point(size=3)

#MDS plot for poisson distance 
mdsPoisData <- data.frame(cmdscale(samplePoisDistMatrix))
mdsPois <- cbind(mdsPoisData, as.data.frame(colData(dds)))
ggplot(mdsPois, aes(X1, X2, color=dex, shape=cell)) +
  geom_point(size=3)


##Differential expression analysis 
dds <- DESeq(dds)

#building results table
(res <- results(dds))

mcols(res, use.names = TRUE)

summary(res)

#lower the false discovery rate threshold 
res.05 <- results(dds, alpha = .05)
table(res.05$padj < 0.5)

#raise the log2 fold change threshold from 0 
  #note: by doing this we test for genes that show significant effects of treatment 
  #on gene counts more than doubling or less than halving 
resLFC1 <- results(dds, lfcThreshold = 1)
table(resLFC1$padj < 0.1)


##other comparisons 
#extract results for the log2 of the fold change of one cell line over another 
results(dds, contrast = c("cell", "N061011", "N61311"))


##multiple testing 
sum(res$pvalue < 0.05, na.rm = TRUE) #5677
sum(!is.na(res$pvalue)) #29391

#genes with adjusted p-value below 10%
sum(res$padj < 0.1, na.rm=TRUE) #4825

#subset the results table to these genes & sort by log2 fold changge to get significant genes w strongest down-regulation
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])

#subset those with strongest up-regulation
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])


##Plotting results
#plotcount
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene=topGene, intgroup=c("dex"))

#custom plots using ggplot function 
geneCounts <- plotCounts(dds, gene=topGene, intgroup=c("dex", "cell"), returnData=TRUE)
ggplot(geneCounts, aes(x=dex, y=count, color=cell)) +
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0), size=3)

ggplot(geneCounts, aes(x=dex, y=count, fill=dex)) +
  scale_y_log10() +
  geom_dotplot(binaxis="y", stackdir="center")

#since deseq2 takes into account the cell line effect the rightmost figure more closely depicts the difference being tested 
ggplot(geneCounts, aes(x=dex, y=count, color=cell, group=cell)) + 
  scale_y_log10() +
  geom_point(size = 3) + geom_line()

#MAplot
plotMA(res, ylim=c(-5,5))

#label individual points on the MAplot
plotMA(resLFC1, ylim=c(-5,5))
topGene <- rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

#histogram of p values 
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50",
     border = "white")

#Gene clustering 
#cluster a subset of the 20 genes with the highest variance across samples 
library("genefilter")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20)

mat <- assay(rld)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("cell", "dex")])
pheatmap(mat, annotation_col = df)
    #Note that a set of genes at the top of the heatmap are separating the N061011 cell line from the others. 
    #In the center of the heatmap, we see a set of genes for which the dex treated samples have higher gene expression

##Independet filtering
#examine ratio of small pvalues ( <0.05) for genes binned by mean normalized count
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean >0], 0:6/6)) #create bins
bins <- cut(resLFC1$baseMean, qs) #bin genes by base mean
levels(bins) <- paste0("-", round(signif(.5*qs[-1] + .5*qs[-length(qs)],2))) #ranme levels of the bins using the middle point 
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < 0.05, na.rm=TRUE)) #calculate the ratio of pvalues <0.05 for each bin 
barplot(ratios, xlab = "mean normalized count", ylab = "ratio of small p values")
  #doing this is useful bc these genes have an influence on the multiple testing adjustment, whoe performance improves if such genes are removed
  #By removing the low count genes from the input to the FDR procedure, we can find more genes to be significant among those that we keep, and so improved the power of our test. This approach is known as independent filtering.

##Annotating and exporting results
library("AnnotationDbi") #helps with mapping various ID schemes to each other 
library("org.Hs.eg.db")

columns(org.Hs.eg.db) #get list of available key types

#add individual columns to results table
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

#now results have the desired external gene IDs 
resOrdered <- res[order(res$padj),]
head(resOrdered)

#exporting results 
resOrderedDF <- as.data.frame(resOrdered)[1:100,]
write.csv(resOrderedDF, file="results.csv")

#alternative method for exporting results
library("ReportingTools") #will automatically generate dynamic HTML documents, 
                          #including links to external databases using gene identifiers and boxplots
                          #summarizing the normalized counts across groups 
htmlRep <- HTMLReport(shortName = "report", title = "My report", reportDirectory = "./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)


##plotting fold changes in genomic space  
#ask for granges output
(resGR <- results(dds, lfcThreshold = 1, format = "GRanges"))
resGR$symbol <- mapIds(org.Hs.eg.db, names(resGR), "SYMBOL", "ENSEMBL")

library("Gviz") #used for plotting the GRanges and associated metadata

#specify a window of 1M bp upstream and downstream from the gene w smallest p-value 
window <- resGR[topGene] + 1e6 
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

#create a vector specifying if the genes in this subset had a low false discovery rate 
sig <- factor(ifelse(resGRsub$padj < 0.1 & !is.na(resGRsub$padj),"sig","notsig"))

#create an axis track specifying our location in the genome
options(ucscChromosomeNames=FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name="gene ranges", feature=sig)
d <- DataTrack(resGRsub, data="log2FoldChange", baseline=0,
                     type="h", name="log2 fold change", strand="+")
#ran dev.off() b/c plot was not working 
plotTracks(list(g,d,a), groupAnnotation="group", notsig="grey", sig="hotpink")



##removing hidden batch effects
#supposing we didn't know there were different cell lines involved in experiment 
# cell line effect on counts would represent hidden and unwanted variation that may be affecting many/all genes in dataset 
library("sva") #statistical methods from this packaged are used to detect such groupings of samples 

#obtain a matrix of normalized counts for which avg count across samples is larger than 1
dat <- counts(dds, normalized=TRUE)
idx <- rowMeans(dat) > 1
dat <- dat[idx,]
mod <- model.matrix(~ dex, colData(dds))
mod0 <- model.matrix(~ 1, colData(dds))
svseq <- svaseq(dat, mod, mod0, n.sv=2) #specify we want to estimate 2 surrogate variables 

svseq$sv

#to see how the sva procedure is albe to identify sources of variation which are correlated with cell line 
par(mfrow=c(2,1),mar=c(3,5,3,1))
stripchart(svseq$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
abline(h=0)
stripchart(svseq$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
abline(h=0)

#add the two surrogate variables as columns to the DESeqDataSet + add them to the design
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + dex

#run the pipeline
ddssva <- DESeq(ddssva)









