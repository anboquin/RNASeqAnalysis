packages <- c("ggplot2", "dplyr")
install.packages(packages, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install ("DESeq2")
BiocManager::install("SummarizedExperiment")
BiocManager::install("pasilla")

library("DESeq2")
library("ggplot2")
library("pasilla")
library("ashr")
library("vsn")

#load in data
pasCts <- system.file("extdata",
                      "pasilla_gene_counts.tsv",
                      package="pasilla", mustWork=TRUE)
pasAnnotation <- system.file("extdata",
                             "pasilla_sample_annotation.csv",
                             package="pasilla", mustWork=TRUE)
cts <- as.matrix(read.csv(pasCts,sep="\t",row.names="gene_id"))
coldata <- read.csv(pasAnnotation, row.names=1)
coldata <- coldata[,c("condition", "type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)

#examine count matrix and column data to see if they're consistent
head(cts, 2)
coldata

#rearrange to remain consistent in terms of sample order
rownames(coldata) <- sub("fb", "", rownames(coldata)) #remove "fb" from row names
all(rownames(coldata) %in% colnames(cts))
##TRUE
all(rownames(coldata) == colnames(cts))
##FALSE
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts)) ##TRUE

##use DESeqDataSetFromMatrix to create dds object
dds_pas <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = ~condition)

#pre-filtering - removing rows w very few reads
keep <- rowSums(counts(dds_pas)) >= 10
dds_pas <- dds_pas[keep,]
dds_pas

#re-leveling
dds_pas$condition <- relevel(dds_pas$condition, ref = "untreated")

#remove levels which do not have samples in current dds
dds_pas$condition <- droplevels(dds_pas$condition)
dds_pas

##Differential expression analysis
dds_pas <- DESeq(dds_pas)
res <- results(dds_pas)
res

#specify coefficient or contrast we want to build a result table for:
res <- results(dds_pas, name="condition_treated_vs_untreated")
res <- results(dds_pas, contrast=c("condition","treated","untreated"))

#more info on results column 
mcols(res)$description

#log fold change shrinkage for visualization and ranking
resultsNames(dds_pas) #intercept condition_treated_vs_untreated
resLFC <- lfcShrink(dds_pas, coef="condition_treated_vs_untreated", type="apeglm")
resLFC

#reorder res table by the smallest p-value
resOrdered <- res[order(res$pvalue),]

summary(res)

#how many adjusted p_values were less than 0.1
sum(res$padj < 0.1, na.rm = TRUE)

#customizing results table
res05 <- results(dds_pas, alpha=0.05)
summary(res05)

#exporting results  to CSV files
write.csv(as.data.frame(resOrdered),
          file="condition_treated_results.csv")
#exporting only results which pass an adujested pbalue above 0.1
resSigPas <- subset(resOrdered, padj < 0.1)
resSigPas


##rerunning analysis using a multi-factor design
#create a copy of the dds
ddsMF <- dds_pas

#change the levels of type to contian only letters
levels(ddsMF$type) #paired-end & single-read
levels(ddsMF$type) <- sub("-.*", "", levels(ddsMF$type))
levels(ddsMF$type) #pared & single

#rerun DESeq
design(ddsMF) <- formula(~ type + condition) #condition is the variable of interest
ddsMF <- DESeq(ddsMF)

#results
resMF <- results(ddsMF)
head(resMF)

resMFType <- results(ddsMF,
                     contrast=c("type", "single", "paired"))
head(resMFType)

#exporting results  to CSV files
write.csv(as.data.frame(resMFType),
          file="MF_condition_treated_results.csv")


#exploring results 
#MA plot
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

#interactive plot 
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

#alternative shrinkage estimators 
resultsNames(dds_pas)

#because we are interested in treated vs untreated we set coef=2
resNorm <- lfcShrink(dds_pas, coef=2, type="normal") 
resAsh <- lfcShrink(dds_pas, coef=2, type="ashr")

#MAplots for shrinkage estimators 
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")

#plot counts - examine counts of reads for a single gene across groups
plotCounts(dds_pas, gene=which.min(res$padj), intgroup="condition")
#customized plot 
d <- plotCounts(dds_pas, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#extracting transformed values 
vsd <- vst(dds_pas, blind=FALSE) #variance stabilizin transformation
rld <- rlog(dds_pas, blind=FALSE) #regularized log transformation 

#effects on transformations on the variance 
# this gives log2(n + 1)
ntd <- normTransform(dds_pas)
meanSdPlot(assay(ntd)) #elevated sd on lower count range
meanSdPlot(assay(vsd))#sd is roughly contant along the whole dynamic range 
meanSdPlot(assay(rld))

#Data quality assessment by sample clustering and visualization 
#heatmap of the count matrix 
library("pheatmap")
select <- order(rowMeans(counts(dds_pas,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_pas)[,c("condition","type")])

#ntd heatmap
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#vsd heatmap
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#rld heatmap 
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap of sample-to-sample distances 
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#principal component plot of samples 
plotPCA(vsd, intgroup=c("condition", "type"))

#customized pca plot
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
