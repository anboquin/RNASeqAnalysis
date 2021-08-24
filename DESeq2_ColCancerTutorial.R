#Griffith Lab tutorial 

#dataset: GXA E-GEOD-50760 
  #data consists of 54 samples from 18 individuals 
  #each has: primary colorectal cancer sample, metastatic liver sample, &  normal sample of surrounding colonic epithilium.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install ("DESeq2")

library("DESeq2")
library("ggplot2")
library("scales") # needed for oob parameter
library("viridis")

#download data 
#read in raw read counts 
rawcounts <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-raw-counts.tsv")
head(rawcounts)

#read in sample mappings
sampledata <- read.delim("http://genomedata.org/gen-viz-workshop/intro_to_deseq2/tutorial/E-GEOD-50760-experiment-design.tsv")
head(sampledata)

#copy
sampledataV2 <- sampledata

#convert count data to a matrix 
geneID <- rawcounts$Gene.ID
sampleIndex <- grepl("SRR\\d+", colnames(rawcounts))
rawcounts <- as.matrix(rawcounts[,sampleIndex])
rownames(rawcounts) <- geneID 
head(rawcounts)

#convert sample variable mappings
head(sampledata)
rownames(sampledata) <- sampledata$Run
keep <- c("Sample.Characteristic.biopsy.site.", "Sample.Characteristic.individual.")
sampledata <- sampledata[,keep]
colnames(sampledata) <- c("tissueType", "individualID")
sampledata$individualID <- factor(sampledata$individualID)
head(sampledata)

#put columns o count data in same order as row names in sample mapping 
rawcounts <- rawcounts[,unique(rownames(sampledata))]
all(colnames(rawcounts) == rownames(sampledata)) #confirm they match 

#rename tissue types 
rename_tissues <- function(x){
  x <- switch(as.character(x), "normal"="normal-looking surrounding colonic epithelium", "primary tumor"="primary colorectal cancer",  "colorectal cancer metastatic in the liver"="metastatic colorectal cancer to the liver")
  return(x)
}
sampledata$tissueType <- unlist(lapply(sampledata$tissueType, rename_tissues))

#reorder: 1) control 2) primary 3) metastatic 
sampledata$tissueType <- factor(sampledata$tissueType, levels = c("normal-looking surrounding colonic epithelium", "primary colorectal cancer", "metastatic colorectal cancer to the liver"))

#create DESeq2Dataset object
dds_cc <- DESeqDataSetFromMatrix(countData=rawcounts, colData=sampledata, design= ~ individualID + tissueType) #individual id = blocking factor & tissue id = comparison variable 

#prefiltering 
dim(dds_cc)
dim(dds_cc[rowSums(counts(dds_cc)) > 5, ]) #check what effect this filter will have 

#apply filter 
dds_cc <- dds_cc[rowSums(counts(dds_cc)) >5,]

#set up multicores
BiocManager::install("BiocParallel")
library(BiocParallel)
register(MulticoreParam(4)) #register number of cores to use 

#running the pipeline
dds_cc <- DESeq(dds_cc)

#extracting restuls 
#tissuetype = perform primary vs normal comparison 
res_cc <- results(dds_cc, contrast=c("tissueType", "primary colorectal cancer", "normal-looking surrounding colonic epithelium"))
summary(res_cc)


#diagnostic plots 
plotMA(res_cc)


# Coerce to a data frame
dds_ResDF <- as.data.frame(res_cc)

# Examine this data frame
head(dds_ResDF)

# Set a boolean column for significance
dds_ResDF$significant <- ifelse(dds_ResDF$padj < .1, "Significant", NA)

# Plot the results similar to DEseq2
ggplot(dds_ResDF, aes(baseMean, log2FoldChange, colour=significant)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + 
  scale_x_log10() + geom_hline(yintercept = 0, colour="tomato1", size=2) + 
  labs(x="mean of normalized counts", y="log fold change") + scale_colour_manual(name="q-value", values=("Significant"="red"), na.value="grey50") + 
  theme_bw()

# add some more detail
  #added overlay density contours 
ggplot(dds_ResDF, aes(baseMean, log2FoldChange, colour=padj)) + 
  geom_point(size=1) + scale_y_continuous(limits=c(-3, 3), oob=squish) + 
  scale_x_log10() + geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + 
  labs(x="mean of normalized counts", y="log fold change") + scale_colour_viridis(direction=-1, trans='sqrt') + 
  theme_bw() + 
  geom_density_2d(colour="black", size=2)
  # we can see that as the average counts increase there is more power to call 
  # a gene as differentially expressed based on the fold change


#viewing normalized counts for a single geneID
#extract counts for gene otop2
otop2Counts <- plotCounts(dds_cc, gene="ENSG00000183034", intgroup=c("tissueType", "individualID"), returnData=TRUE)

#plot the data 
colourPallette <- c("#7145cd","#bbcfc4","#90de4a","#cd46c1","#77dd8e","#592b79","#d7c847","#6378c9","#619a3c","#d44473","#63cfb6","#dd5d36","#5db2ce","#8d3b28","#b1a4cb","#af8439","#c679c0","#4e703f","#753148","#cac88e","#352b48","#cd8d88","#463d25","#556f73")
ggplot(otop2Counts, aes(x=tissueType, y=count, colour=individualID, group=individualID)) + 
  geom_point() + 
  geom_line() + 
  theme_bw() + 
  theme(axis.text.x=element_text(angle=15, hjust=1)) + 
  scale_colour_manual(values=colourPallette) + 
  guides(colour=guide_legend(ncol=3)) + ggtitle("OTOP2")
    #almost all individuals show down-regulation of this gene in both 
    #priamry tumor and metastasis samples compared to the normal 

#visualizing expression data with a heatmap 
#Transform count data using the variance stablilizing transform
dds_2VST <- vst(dds_cc)

#convert the DESeq transformed object to a data frame
dds_2VST <- assay(dds_2VST) #extract transformed vaues 
dds_2VST <- as.data.frame(dds_2VST)
dds_2VST$Gene <- rownames(dds_2VST)
head(dds_2VST)

#keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(dds_ResDF[dds_ResDF$padj <= .05 & abs(dds_ResDF$log2FoldChange) > 3,])
dds_2VST <- dds_2VST[dds_2VST$Gene %in% sigGenes,]

#convert the VST counts to long format for ggplot2
library(reshape2)

#compare wide vs long version
dds2VST_wide <- dds_2VST
dds2VST_long <- melt(dds_2VST, id.vars=c("Gene")) #long is required for ggplot2

head(dds2VST_wide)
head(dds2VST_long)

#overwrite  original data frame with the long format
dds_2VST <- melt(dds_2VST, id.vars=c("Gene"))

#Make a heatmap
heatmap <- ggplot(dds_2VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + 
  scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap


#clustering 
#convert the significant genes back to a matrix for clustering
dds2VSTMatrix <- dcast(dds_2VST, Gene ~ variable) #convert transfromed values back to wide format 
rownames(dds2VSTMatrix) <- dds2VSTMatrix$Gene
dds2VSTMatrix$Gene <- NULL

#compute a distance calculation on both dimensions of the matrix
distanceGene <- dist(dds2VSTMatrix) 
distanceSample <- dist(t(dds2VSTMatrix)) #t = transpose 

#cluster based on the distance calculations
clusterGene <- hclust(distanceGene, method="average") #hierarchical clustering
clusterSample <- hclust(distanceSample, method="average")

#construct a dendogram for samples
install.packages("ggdendro")
library(ggdendro)
sampleModel <- as.dendrogram(clusterSample)
sampleDendrogramData <- segment(dendro_data(sampleModel, type = "rectangle"))
sampleDendrogram <- ggplot(sampleDendrogramData) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  theme_dendro()

#re-factor samples for ggplot2
dds_2VST$variable <- factor(dds_2VST$variable, levels=clusterSample$labels[clusterSample$order])

#construct the heatmap. note that at this point we have only clustered the samples NOT the genes
heatmap <- ggplot(dds_2VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + 
  scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())
heatmap

#combine the dendrogram and the heatmap
install.packages("gridExtra")
library(gridExtra)
grid.arrange(sampleDendrogram, heatmap, ncol=1, heights=c(1,5))
    #plot widths dont match up b/c we have a legend in one plot but not in the other 

##fixing the plot
library(gtable) #allows us to view and manipulate grobs as tables 
library(grid)

#modify ggplot objects
sampleDendrogram_1 <- sampleDendrogram + 
  scale_x_continuous(expand=c(.0085, .0085)) + 
  scale_y_continuous(expand=c(0, 0))
heatmap_1 <- heatmap + 
  scale_x_discrete(expand=c(0, 0)) + 
  scale_y_discrete(expand=c(0, 0))

#convert both grid based objects to grobs to make them easier to manipulate 
sampleDendrogramGrob <- ggplotGrob(sampleDendrogram_1)
heatmapGrob <- ggplotGrob(heatmap_1)

#check witdths of each grob
sampleDendrogramGrob$widths
heatmapGrob$widths

#add in missing columns 
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[7], 6)
sampleDendrogramGrob <- gtable_add_cols(sampleDendrogramGrob, heatmapGrob$widths[8], 7)

#make sure every width between the two grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths) #find max widths for each grob
sampleDendrogramGrob$widths <- as.list(maxWidth) #override each grobs width with max widths 
heatmapGrob$widths <- as.list(maxWidth)

#arrange grobs into a plot 
finalGrob <- arrangeGrob(sampleDendrogramGrob, heatmapGrob, ncol=1, heights=c(2,5))

#draw the plot
grid.draw(finalGrob)

#add plot b/w dendrogram and heatmap to show the tissue type 
#re-order sample data to match clustering 
sampledataV2$Run <- factor(sampledataV2$Run, levels=clusterSample$labels[clusterSample$order])

#construct a plot to show the clinical data
colours <- c("#743B8B", "#8B743B", "#8B3B52")
sampleClinical <- ggplot(sampledataV2, aes(x=Run, y=1, fill=Sample.Characteristic.biopsy.site.)) + geom_tile() + scale_x_discrete(expand=c(0, 0)) + scale_y_discrete(expand=c(0, 0)) + scale_fill_manual(name="Tissue", values=colours) + theme_void()

#convert the clinical plot to a grob
sampleClinicalGrob <- ggplotGrob(sampleClinical)

#make sure every width between all grobs is the same
maxWidth <- unit.pmax(sampleDendrogramGrob$widths, heatmapGrob$widths, sampleClinicalGrob$widths)
sampleDendrogramGrob$widths <- as.list(maxWidth)
heatmapGrob$widths <- as.list(maxWidth)
sampleClinicalGrob$widths <- as.list(maxWidth)

#arrange and output the final plot
finalGrob <- arrangeGrob(sampleDendrogramGrob, sampleClinicalGrob, heatmapGrob, ncol=1, heights=c(2,1,5))
grid.draw(finalGrob)




