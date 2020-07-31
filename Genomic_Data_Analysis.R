SRR6809855 <- read.delim("SRR6809855_GRCh38.tsv")
SRR6809856 <- read.delim("SRR6809856_GRCh38.tsv")
SRR6809862 <- read.delim("SRR6809862_GRCh38.tsv")
SRR6809883 <- read.delim("SRR6809883_GRCh38.tsv")
SRR6809886 <- read.delim("SRR6809886_GRCh38.tsv")
SRR6809890 <- read.delim("SRR6809890_GRCh38.tsv")


plot(SRR6809855$tpm, SRR6809883$tpm)


head(SRR6809855)
head(subset(SRR6809855, est_counts > 0))

ped <- data.frame(
  pedTB_samp1 = round(SRR6809855$est_counts),
  pedTB_samp2 = round(SRR6809856$est_counts),
  pedTB_samp3 = round(SRR6809862$est_counts),
  HC_samp1 = round(SRR6809883$est_count),
  HC_samp2 = round(SRR6809886$est_count),
  HC_samp3 = round(SRR6809890$est_count)
)

rownames(ped) <- SRR6809855$target_id

head(ped)
head(subset(ped, pedTB_samp1>0))



coldata <- data.frame(
  treatment = c(rep("pedTB",3), rep("HC",3))
)

rownames(coldata) <- colnames(ped)
coldata


library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = ped,
                              colData = coldata,
                              design = ~treatment)
dds <- DESeq(dds)
res <- results(dds)
res_full <- results(dds)
summary(res_full)

#perform minimum pre-filtering 
dds_filt <- dds[rowSums(counts(dds)) >= 10,]
dds_filtered <- DESeq(dds_filt)
res_filt <- results(dds_filtered)
summary(res_filt)
summary(res)
subset(res, baseMean>0)
res <- na.omit(res)
summary(res)
res
res_filt2 <- na.omit(res_filt)
summary(res_filt2)
#remove NAs

#look at gene with highest LFC
res_ordered <- res_filt2[order(-res_filt2$log2FoldChange),]
high_lfc <- res_ordered[1,]

#look at gene with lowest p-value
low_p <- res_filt2[which.min(res_filt2$padj),]


#look at summary for highlfc gene and lowp gene
summary(high_lfc)
summary(low_p)

plotCounts(dds_filtered, gene = "ENST00000424832.6", intgroup = "treatment")

plotCounts(dds_filtered, gene = "ENST00000452392.2", intgroup = "treatment")


library(apeglm)
resultsNames(dds_filtered)
resLFC <- lfcShrink(dds_filtered,
                    coef="treatment_pedTB_vs_HC",
                    type="apeglm")
resLFC <- na.omit(resLFC)
plot(res_filt2$log2FoldChange, resLFC$log2FoldChange,
     xlab="unshrunk LFC", ylab="shrunk LFC")

plot(resLFC$baseMean, resLFC$log2FoldChange)
#res[order(res$padj),]
#res[order(res$log2FoldChange),]
plotMA(res_filt2, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")



dev.off()
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment", returnData=TRUE) 
library("ggplot2") 
ggplot(d, aes(x=treatment, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(1,5,10))


ntd <- normTransform(dds_filtered)
vsd <- vst(dds_filtered, blind=FALSE) 
rld <- rlog(dds_filtered, blind=FALSE)
head(assay(vsd), 3)

colData(dds_filtered)

dds_filtered$treatment
library("pheatmap") 
select <- order(rowMeans(counts(dds_filtered,normalized=TRUE)), decreasing=TRUE)[1:20] 
annot <- dds_filtered$treatment
pheatmap(assay(vsd)[select,],
         cluster_rows = FALSE,
         show_rownames=FALSE,
         cluster_cols=FALSE,
         annotation_col=df)
annotations <- c(rep("pedTB",3), rep("HC",3))
df <- as.data.frame(colData(dds_filtered)[,c("treatment")])
pheatmap(assay(vsd)[select,],cluster_cols=FALSE,labels_col=annotations)

pheatmap(assay(vsd)[select,])
#get order of genes as seen on heatmap
h_m<- pheatmap(assay(vsd)[select,])
hc <- as.hclust( h_m$tree_row)
cutree( hc, h=11 )[hc$order] #tb 
cutree( hc, h=10 )[hc$order] #hc
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists) 


rownames(sampleDistMatrix) <- paste(vsd$treatment) 


colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(7,"Blues")))(255) 
colors <- colorRampPalette(c("blue","black","red"))(n=600)

pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists)




library("RColorBrewer")



plotPCA(vsd,intgroup=("treatment"))


