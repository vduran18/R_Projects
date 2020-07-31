library(readr)
library(ggplot2)

SRR493372 <- read_tsv("/Users/valeriaduran/Library/Mobile Documents/com~apple~CloudDocs/MATH6388/abundance_SRR493372.tsv")
SRR493374 <- read_tsv("/Users/valeriaduran/Library/Mobile Documents/com~apple~CloudDocs/MATH6388/abundance_SRR493374.tsv")

head(SRR493372)

#merge data 

kallipso <- merge(SRR493372, SRR493374, by = 'target_id')

ggplot(kallipso, aes(tpm.x, tpm.y)) +geom_point() +theme_minimal() 
  
ggplot(kallipso, aes(x=tpm.x, y=tpm.y, color = target_id)) + 
         geom_point()
dev.off()
ggplot() 
  
ggplot(kallipso, aes(target_id, y = value, color = variable)) + 
  geom_point(aes(y = tpm.x, col = "tpm.x")) + 
  geom_point(aes(y = tpm.y, col = "tpm.y")) + 
  labs(xlab = 




plot(SRR493372$tpm, SRR493374$tpm)




BiocManager::install("airway")
library(airway) 
data(airway)
library("DESeq2")
ddsSE <- DESeqDataSet(airway, design = ~ cell + dex)

colData(ddsSE)
ddsSE <- DESeq(ddsSE)
resSE <- results(ddsSE)
dds <- ddsSE[rowSums(counts(ddsSE)) >= 10,]
dds <- DESeq(dds)

res <- results(dds)

summary(resSE)
summary(res)

res_ordered <- res[order(-res$log2FoldChange),]
res_ordered[1,]
rm(random)
summary(res[which.min(res$padj),])

random <- res[sample(nrow(res),1),]
random
summary(random)
res[10,]


plotCounts(dds, gene = "ENSG00000243444", intgroup = "dex")

DGEgenes <- rownames(random)
resSE[DGEgenes,]
summary(random)
plotMA(res)




library(apeglm)
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="dex_untrt_vs_trt",
                    type="apeglm")
resLFC



ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("cell","dex")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE,
         show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)



plotPCA(vsd, intgroup=c("cell", "dex"))



library("RColorBrewer")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists) 
rownames(sampleDistMatrix) <- paste(vsd$cell, vsd$dex, sep="-") 
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) 
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)




