packages <- c("ggplot2", "dplyr")
install.packages(packages, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install ("DESeq2")
BiocManager::install("SummarizedExperiment")
BiocManager::install("ashr")

library("DESeq2")
library("ggplot2")
library("ashr")
library("vsn")

#load data 
airCounts <- read.csv("airway_scaledcounts.csv", header = TRUE, sep=",",row.names = "ensgene")
dim(airCounts)
head(airCounts)  #gene identifies + expression values for each gene
airCounts <- as.matrix(airCounts) #convert to matrix

airmetadata <- read.csv("airway_metadata.csv",header = TRUE, sep=",", row.names = 1)
airmetadata <- airmetadata[,c("dex", "celltype")] #removed geo_id column
head(airmetadata)
airmetadata$dex <- factor(airmetadata$dex)
airmetadata$celltype <- factor(airmetadata$celltype)

#check order with respect to samples
all(rownames(airmetadata) %in% colnames(airCounts)) #TRUE
all(rownames(airmetadata) == colnames(airCounts)) #true 


#build deseqdataset object
dds_air <-DESeqDataSetFromMatrix(countData = airCounts,
                                colData = airmetadata,
                                design = ~dex)
dds_air #dim 38694 8 
head(assay(dds_air))

#pre-filtering
keep <- rowSums(counts(dds_air)) >=10
dds_air <- dds_air[keep,]
dds_air # dim: 19271 8 

#releveling 
dds_air$dex <- relevel(dds_air$dex, ref="control")

#running pipeline
dds_air <- DESeq(dds_air)
dds_air

#results table
results_air <- results(dds_air)
results_air
summary(results_air)

#sort by smallest p-value
res_airOrdered <- results_air[order(results_air$pvalue),]
head(res_airOrdered)

#how many adjusted p-values less than 0.1?
sum(results_air$padj < 0.1, na.rm = TRUE) 
  #4041

#how many adjusted p-values less than 0.05?
res_air05 <- results(dds_air, alpha=0.05)
summary(res_air05)
sum(res_air05$padj <0.05, na.rm = TRUE)
  #3324


#exporting results
write.csv(as.data.frame(res_airOrdered),
          file = "airwaytreated_results.csv")
#exporting only results which pass an adjusted pvalue above 0.1
resSig <- subset(res_airOrdered, padj < 0.1)

#log fold change shrinkage
resLFC_air <- lfcShrink(dds_air, coef="dex_treated_vs_control", type = "apeglm")
resLFC_air

#alternative shrinkage estimators
resairNorm <- lfcShrink(dds_air, coef = 2, type = "normal")
resairAsh <- lfcShrink(dds_air, coef = 2, type = "ashr")

#Exploring results 
#MA plots
plotMA(results_air, ylim=c(-3,3))
plotMA(resLFC_air, ylim=c(-3,3))

#MAplots for schrinkage estimators 
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC_air, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resairNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resairAsh, xlim=xlim, ylim=ylim, main="ashr")


#plot counts 
#examine plot for gene with samllest adjusted pvalue
plotCounts(dds_air, gene=which.min(results_air$padj), intgroup="dex")
#customized plot 
d_air <- plotCounts(dds_air, gene = which.min(results_air$padj), intgroup = "dex", returnData = TRUE)
ggplot(d_air, aes(x=dex, y = count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))

#examine plots bw treated and control groups for top 6 genes 
par(mfrow=c(2,3))
plotCounts(dds_air, gene="ENSG00000152583", intgroup="dex")
plotCounts(dds_air, gene="ENSG00000179094", intgroup="dex")
plotCounts(dds_air, gene="ENSG00000116584", intgroup="dex")
plotCounts(dds_air, gene="ENSG00000189221", intgroup="dex")
plotCounts(dds_air, gene="ENSG00000120129", intgroup="dex")
plotCounts(dds_air, gene="ENSG00000148175", intgroup="dex")


#volcano plot 
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(results_air, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

#Add coloured points: blue if padj <0.1, red if log2FC>1 and padj<0.05
with(subset(results_air, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(results_air, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


##Extracting ransformed values 
#vst transformation 
vsdata <- vst(dds_air, blind = FALSE) 
plotPCA(vsdata, intgroup="dex") #principal component analysis 
                                #look at how samples group by treatment
plotPCA(vsdata, intgroup=c("dex", "celltype"))

#customized PCA plot for vsdata
pcaData <- plotPCA(vsdata, intgroup=c("dex", "celltype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dex, shape=celltype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

#rlog transformation
rldata <- rlog(dds_air, blind = FALSE)
plotPCA(rldata, intgroup="dex")
plotPCA(rldata, intgroup=c("dex", "celltype"))

#customized PCA plot for rldata
pcaData <- plotPCA(rldata, intgroup=c("dex", "celltype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=dex, shape=celltype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()

#plot sd of transformed data across samples against mean using the shifted log transformation
ntd <- normTransform(dds_air) #this gives log2(n+1)
meanSdPlot(assay(ntd))

#plot using the vst
meanSdPlot(assay(vsdata))

#plot using rlog transformation
meanSdPlot(assay(rldata))

#constructing heatmaps 
library("pheatmap")
select <- order(rowMeans(counts(dds_air,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds_air)[,c("dex","celltype")])

#heatmap using shifted log transformation
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap using vst
pheatmap(assay(vsdata)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap using rlog
pheatmap(assay(rldata)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

#heatmap of sample-to-sample distances  - sample clustering
#apply dist function to the transpose of the transformed count matrix to get sample-to-sample distances
sampleDists <- dist(t(assay(vsdata)))


library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdata$dex, vsdata$celltype)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


